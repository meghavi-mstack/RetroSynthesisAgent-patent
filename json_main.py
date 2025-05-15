# --- Add this block at the TOP of json_main.py ---
import ssl
import os
import json
import re
import traceback
import argparse
import fitz  # PyMuPDF
import pubchempy

# WARNING: Disables SSL certificate verification globally!
# Use only if necessary and you understand the risks (e.g., trusted network).
if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
        getattr(ssl, '_create_unverified_context', None)):
    print("!!! WARNING: SSL CERTIFICATE VERIFICATION DISABLED GLOBALLY !!!")
    ssl._create_default_https_context = ssl._create_unverified_context
# --- End of SSL bypass block ---

from RetroSynAgent.treeBuilder import Tree, TreeLoader
from RetroSynAgent.pdfProcessor import PDFProcessor
from RetroSynAgent.knowledgeGraph import KnowledgeGraph
from RetroSynAgent import prompts
from RetroSynAgent.GPTAPI import GPTAPI
from RetroSynAgent.patentPDFDownloader import PatentPDFDownloader
from RetroSynAgent.pdfDownloader import PDFDownloader
from RetroSynAgent.name_to_smiles_fixed import NameToSMILES
from RetroSynAgent.entityAlignment import EntityAlignment
from RetroSynAgent.treeExpansion import TreeExpansion
from RetroSynAgent.reactionsFiltration import ReactionsFiltration


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Process PDFs and extract reactions.")
    parser.add_argument('--material', type=str, required=True,
                        help="Chemical name or SMILES string of the target molecule (e.g., 'aspirin' or 'CC(=O)OC1=CC=CC=C1C(=O)O').")
    parser.add_argument('--num_results', type=int, required=True,
                        help="Maximum number of patent PDFs to download.")
    parser.add_argument('--alignment', type=str, default="False", choices=["True", "False"],
                        help="Whether to align entities except for root node.")
    parser.add_argument('--expansion', type=str, default="False", choices=["True", "False"],
                        help="Whether to expand the tree with additional literature.")
    parser.add_argument('--filtration', type=str, default="False", choices=["True", "False"],
                        help="Whether to filter reactions.")
    parser.add_argument('--retrieval_mode', type=str, default="patent-patent",
                        choices=["patent-patent", "paper-paper", "both-both"],
                        help="Document retrieval mode: patent-patent (patents for both), paper-paper (papers for both), both-both (both patents and papers for both initial and expansion)")
    parser.add_argument('--output', type=str, default=None,
                        help="Output JSON file path (defaults to [material]_pathways.json if not specified)")
    return parser.parse_args()


def countNodes(tree):
    node_count = tree.get_node_count()
    return node_count


def searchPathways(tree):
    all_path = tree.find_all_paths()
    return all_path


def main(material,
         num_results,
         alignment,
         expansion,
         filtration,
         retrieval_mode="patent-paper",
         output_file=None):
    # If output_file is not specified, use the material name
    if output_file is None:
        output_file = f"{material.lower()}_pathways.json"
    try:
        print("Starting json_main function...")
        print(f"Material: {material}")
        print(f"Number of results: {num_results}")
        print(f"Alignment: {alignment}")
        print(f"Expansion: {expansion}")
        print(f"Filtration: {filtration}")
        print(f"Retrieval mode: {retrieval_mode}")
        print(f"Output file: {output_file}")

        pdf_folder_name = 'pdf_pi'
        result_folder_name = 'res_pi'
        result_json_name = 'llm_res'
        tree_folder_name = 'tree_pi'
        os.makedirs(tree_folder_name, exist_ok=True)
        os.makedirs(result_folder_name, exist_ok=True)
        entityalignment = EntityAlignment()
        treeloader = TreeLoader()
        tree_expansion = TreeExpansion()
        reactions_filtration = ReactionsFiltration()

        ### extractInfos

        # 1  Convert chemical name to SMILES if needed
        print(f"Processing material: {material}")
        smiles = material
        valid_smiles = True

        # Improved detection of chemical names vs SMILES strings
        # Common chemical name patterns (like brackets in nomenclature)
        chemical_name_patterns = [
            r'\[[0-9\.]+\.[0-9\.]+\.[0-9\.]+\]',  # Matches patterns like [1.1.1] in Bicyclo[1.1.1]pentane
            r'\([a-zA-Z0-9,\-]+\)',  # Matches patterns like (R)-, (S)-, (E)-, (Z)-, etc.
            r'[a-zA-Z]{5,}',  # Long alphabetic strings are likely chemical names
            r'[0-9]+-[a-zA-Z]',  # Numbered prefixes like 1,2-di, 2,3,4-tri
            r'(di|tri|tetra|penta|hexa|hepta|octa|nona|deca)[a-z]'  # Common prefixes
        ]

        # SMILES-specific patterns that are unlikely in chemical names
        smiles_specific_patterns = [
            r'[=#@\\/]',  # SMILES operators
            r'\[\w+\]',  # Atom specifications like [C@H], [NH2], etc.
            r'c1[cnosp]',  # Aromatic carbon patterns
            r'C[0-9]\(',  # Ring specifications
            r'\(=O\)',  # Carbonyl groups in SMILES
            r'\.[A-Z]',  # Disconnected structures
        ]

        # Check if it matches chemical name patterns
        is_likely_chemical_name = any(re.search(pattern, material) for pattern in chemical_name_patterns)

        # Check if it contains SMILES-specific patterns
        has_smiles_patterns = any(re.search(pattern, material) for pattern in smiles_specific_patterns)

        # Additional check for common chemical name words
        common_chemical_words = ['acid', 'amine', 'ether', 'alcohol', 'aldehyde', 'ketone', 'ester',
                                 'benzene', 'methyl', 'ethyl', 'propyl', 'butyl', 'pentyl', 'hexyl',
                                 'cyclo', 'bicyclo', 'tricyclo', 'cubane', 'adamantane', 'fullerene']
        contains_chemical_word = any(word.lower() in material.lower() for word in common_chemical_words)

        # Special case for Bicyclo compounds which are often misidentified
        is_bicyclo_compound = 'bicyclo' in material.lower() and '[' in material and ']' in material

        # Decision logic: prioritize chemical name detection
        if is_likely_chemical_name or contains_chemical_word or is_bicyclo_compound:
            print(f"Input appears to be a chemical name. Converting to SMILES...")
            success, result = NameToSMILES.convert(material)
            if success:
                smiles = result
                print(f"Successfully converted '{material}' to SMILES: {smiles}")
            else:
                # For Bicyclo compounds, try a direct lookup for common structures
                if is_bicyclo_compound:
                    bicyclo_smiles_map = {
                        'bicyclo[1.1.1]pentane': 'C1C2CC1C2',
                        'bicyclo[2.2.1]heptane': 'C1CC2CCC1C2',
                        'bicyclo[2.2.2]octane': 'C1CC2CCC1CC2',
                        'bicyclo[3.3.1]nonane': 'C1CC2CCCC(C1)C2',
                        'bicyclo[4.4.0]decane': 'C1CCC2CCCCC2C1'
                    }

                    # Try case-insensitive match
                    success = False
                    for name, smile in bicyclo_smiles_map.items():
                        if name.lower() == material.lower():
                            smiles = smile
                            print(f"Found SMILES for '{material}' in lookup table: {smiles}")
                            success = True
                            break

                    if not success:
                        print(f"Warning: Could not convert '{material}' to SMILES: {result}")
                        print(f"Proceeding with the original input as SMILES.")
                        valid_smiles = False
                else:
                    print(f"Warning: Could not convert '{material}' to SMILES: {result}")
                    print(f"Proceeding with the original input as SMILES.")
                    valid_smiles = False
        else:
            print(f"Input appears to be a SMILES string already. Proceeding without conversion.")
            # Validate if it's actually a valid SMILES string
            try:
                # Simple validation check - if it contains brackets, make sure they're balanced
                if '[' in material and material.count('[') != material.count(']'):
                    print(f"Warning: Unbalanced brackets in SMILES string. This might not be a valid SMILES.")
                    print(f"Attempting to convert as a chemical name instead...")
                    success, result = NameToSMILES.convert(material)
                    if success:
                        smiles = result
                        print(f"Successfully converted '{material}' to SMILES: {smiles}")
                    else:
                        # Even if conversion fails, we'll still treat it as a SMILES string
                        # since it was detected as one initially
                        valid_smiles = True
                        print(f"Proceeding with the original input as SMILES.")
                else:
                    # If brackets are balanced or there are no brackets, it's likely a valid SMILES
                    valid_smiles = True
            except Exception as e:
                print(f"Warning: Error validating SMILES string: {str(e)}")
                # Even if there's an error, we'll still try to use it as a SMILES
                valid_smiles = True
                print(f"Proceeding with the original input as SMILES despite validation error.")

        # Decide on initial retrieval method based on valid SMILES and retrieval mode
        pdf_name_list = []
        os.makedirs(pdf_folder_name, exist_ok=True)

        if retrieval_mode == "both-both":
            # Retrieve from both patents and papers
            # Calculate how many documents to retrieve from each source
            # Divide the requested number of results between patents and papers
            patents_to_retrieve = num_results // 2
            papers_to_retrieve = num_results - patents_to_retrieve

            print(f"Retrieving {patents_to_retrieve} patents and {papers_to_retrieve} papers for initial retrieval...")

            # First retrieve from patents if we have a valid SMILES
            patent_pdf_list = []
            if valid_smiles:
                try:
                    print("Initializing PatentPDFDownloader for initial retrieval...")
                    downloader = PatentPDFDownloader(pdf_folder_name=pdf_folder_name, max_patents=patents_to_retrieve)
                    print("Calling process_smile method...")
                    # Default Redis connection parameters
                    redis_host = os.getenv("REDIS_HOST", "localhost")
                    redis_port = int(os.getenv("REDIS_PORT", 6379))
                    redis_db = int(os.getenv("REDIS_DB", 0))
                    patent_pdf_list = downloader.process_smile(smiles, redis_host, redis_port, redis_db)
                    print(f'Successfully downloaded {len(patent_pdf_list)} PDFs from patents for SMILES: {smiles}')
                except ValueError as e:
                    print(f"Error with patent search: {str(e)}")
                    print("Will still proceed with academic paper search.")
            else:
                print("No valid SMILES available for patent search. Will only retrieve academic papers.")

            # Then retrieve from academic papers
            print("Initializing PDFDownloader for academic paper retrieval...")
            downloader = PDFDownloader(material, pdf_folder_name=pdf_folder_name, num_results=papers_to_retrieve, n_thread=3)
            paper_pdf_list = downloader.main()
            print(f'Successfully downloaded {len(paper_pdf_list)} PDFs from academic papers for {material}')

            # Combine the results
            pdf_name_list = patent_pdf_list + paper_pdf_list
            print(f'Total PDFs downloaded: {len(pdf_name_list)} ({len(patent_pdf_list)} patents, {len(paper_pdf_list)} papers)')

        elif retrieval_mode.startswith("patent") and valid_smiles:
            # Patent-based initial retrieval (requires valid SMILES)
            try:
                print("Initializing PatentPDFDownloader for initial retrieval...")
                downloader = PatentPDFDownloader(pdf_folder_name=pdf_folder_name, max_patents=num_results)
                print("Calling process_smile method...")
                # Default Redis connection parameters
                redis_host = os.getenv("REDIS_HOST", "localhost")
                redis_port = int(os.getenv("REDIS_PORT", 6379))
                redis_db = int(os.getenv("REDIS_DB", 0))
                pdf_name_list = downloader.process_smile(smiles, redis_host, redis_port, redis_db)
                print(f'Successfully downloaded {len(pdf_name_list)} PDFs from patents for SMILES: {smiles}')
            except ValueError as e:
                print(f"Error with patent search: {str(e)}")
                print("Falling back to academic paper search.")
                # Fall back to academic paper search
                print("Initializing PDFDownloader for initial retrieval (fallback)...")
                downloader = PDFDownloader(material, pdf_folder_name=pdf_folder_name, num_results=num_results, n_thread=3)
                pdf_name_list = downloader.main()
                print(f'Successfully downloaded {len(pdf_name_list)} PDFs from academic papers for {material}')
        else:  # paper-based initial retrieval
            print("Initializing PDFDownloader for initial retrieval...")
            downloader = PDFDownloader(material, pdf_folder_name=pdf_folder_name, num_results=num_results, n_thread=3)
            pdf_name_list = downloader.main()
            print(f'Successfully downloaded {len(pdf_name_list)} PDFs from academic papers for {material}')

        if not pdf_name_list:
            print("No PDFs were downloaded. This could be because:")
            print("1. No patent IDs were found in Redis for this SMILES string")
            print("2. The patents found did not have downloadable PDFs")
            print("3. There was an error connecting to the patent database")
            print("\nPlease check that Redis is properly set up with patent data.")
            print("You can run the debug_downloader.py script to diagnose Redis issues.")
            return {"error": "No PDFs downloaded"}

        # 2 Extract infos from PDF about reactions
        pdf_processor = PDFProcessor(pdf_folder_name=pdf_folder_name, result_folder_name=result_folder_name,
                                     result_json_name=result_json_name)
        pdf_processor.load_existing_results()
        pdf_processor.process_pdfs_txt(save_batch_size=2)

        ### treeBuildWOExapnsion
        results_dict = entityalignment.alignRootNode(result_folder_name, result_json_name, material)

        # 4 construct kg & tree
        tree_name_wo_exp = tree_folder_name + '/' + material + '_wo_exp.pkl'
        if not os.path.exists(tree_name_wo_exp):
            tree_wo_exp = Tree(material.lower(), result_dict=results_dict)
            print('Starting to construct RetroSynthetic Tree...')
            tree_wo_exp.construct_tree()
            treeloader.save_tree(tree_wo_exp, tree_name_wo_exp)
        else:
            tree_wo_exp = treeloader.load_tree(tree_name_wo_exp)
            print('RetroSynthetic Tree wo expansion already loaded.')
        node_count_wo_exp = countNodes(tree_wo_exp)
        all_path_wo_exp = searchPathways(tree_wo_exp)
        print(f'The tree contains {node_count_wo_exp} nodes and {len(all_path_wo_exp)} pathways before expansion.')

        if alignment:
            print('Starting to align the nodes of RetroSynthetic Tree...')

            ### WO Expansion
            tree_name_wo_exp_alg = tree_folder_name + '/' + material + '_wo_exp_alg.pkl'
            if not os.path.exists(tree_name_wo_exp_alg):
                reactions_wo_exp = tree_wo_exp.reactions
                reactions_wo_exp_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_wo_exp)
                reactions_wo_exp_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_wo_exp_alg_1)
                tree_wo_exp_alg = Tree(material.lower(), reactions=reactions_wo_exp_alg_all)
                tree_wo_exp_alg.construct_tree()
                treeloader.save_tree(tree_wo_exp_alg, tree_name_wo_exp_alg)
            else:
                tree_wo_exp_alg = treeloader.load_tree(tree_name_wo_exp_alg)
                print('aligned RetroSynthetic Tree wo expansion already loaded.')
            node_count_wo_exp_alg = countNodes(tree_wo_exp_alg)
            all_path_wo_exp_alg = searchPathways(tree_wo_exp_alg)
            print(
                f'The aligned tree contains {node_count_wo_exp_alg} nodes and {len(all_path_wo_exp_alg)} pathways before expansion.')
            tree_wo_exp = tree_wo_exp_alg  # Update tree_wo_exp for further processing

        ## treeExpansion
        # Skip expansion if not requested
        if expansion:
            # 5 kg & tree expansion
            results_dict_additional = None
            results_dict_additional = tree_expansion.treeExpansion(result_folder_name, result_json_name,
                                                                  results_dict, material, expansion=True, max_iter=5,
                                                                  retrieval_mode=retrieval_mode, smiles=smiles)
            if results_dict_additional:
                results_dict = tree_expansion.update_dict(results_dict, results_dict_additional)
                print(f"Added {len(results_dict_additional)} additional reaction entries from expansion.")

            tree_name_exp = tree_folder_name + '/' + material + '_w_exp.pkl'
            if not os.path.exists(tree_name_exp):
                tree_exp = Tree(material.lower(), result_dict=results_dict)
                print('Starting to construct Expanded RetroSynthetic Tree...')
                tree_exp.construct_tree()
                treeloader.save_tree(tree_exp, tree_name_exp)
            else:
                tree_exp = treeloader.load_tree(tree_name_exp)
                print('RetroSynthetic Tree w expansion already loaded.')
        else:
            # Use the non-expanded tree if expansion is not requested
            tree_exp = tree_wo_exp
            print('Skipping tree expansion as requested.')

        # nodes & pathway count (tree w exp)
        node_count_exp = countNodes(tree_exp)
        all_path_exp = searchPathways(tree_exp)
        print(f'The tree contains {node_count_exp} nodes and {len(all_path_exp)} pathways after expansion.')

        if alignment and expansion:
            ### Expansion alignment (only if both alignment and expansion are enabled)
            tree_name_exp_alg = tree_folder_name + '/' + material + '_w_exp_alg.pkl'
            if not os.path.exists(tree_name_exp_alg):
                reactions_exp = tree_exp.reactions
                reactions_exp_alg_1 = entityalignment.entityAlignment_1(reactions_dict=reactions_exp)
                reactions_exp_alg_all = entityalignment.entityAlignment_2(reactions_dict=reactions_exp_alg_1)
                tree_exp_alg = Tree(material.lower(), reactions=reactions_exp_alg_all)
                tree_exp_alg.construct_tree()
                treeloader.save_tree(tree_exp_alg, tree_name_exp_alg)
            else:
                tree_exp_alg = treeloader.load_tree(tree_name_exp_alg)
                print('aligned RetroSynthetic Tree w expansion already loaded.')
            node_count_exp_alg = countNodes(tree_exp_alg)
            all_path_exp_alg = searchPathways(tree_exp_alg)
            print(
                f'The aligned tree contains {node_count_exp_alg} nodes and {len(all_path_exp_alg)} pathways after expansion.')
            tree_exp = tree_exp_alg  # Update tree_exp for further processing

        all_pathways_w_reactions = reactions_filtration.getFullReactionPathways(tree_exp)

        ## Filtration
        if filtration:
            # filter reactions based on conditions
            reactions_txt_filtered = reactions_filtration.filterReactions(tree_exp)
            # build & save tree
            tree_name_filtered = tree_folder_name + '/' + material + '_filtered' + '.pkl'
            if not os.path.exists(tree_name_filtered):
                print('Starting to construct Filtered RetroSynthetic Tree...')
                tree_filtered = Tree(material.lower(), reactions_txt=reactions_txt_filtered)
                tree_filtered.construct_tree()
                treeloader.save_tree(tree_filtered, tree_name_filtered)
            else:
                tree_filtered = treeloader.load_tree(tree_name_filtered)
                print('Filtered RetroSynthetic Tree already loaded.')
            node_count_filtered = countNodes(tree_filtered)
            all_path_filtered = searchPathways(tree_filtered)
            print(
                f'The tree contains {node_count_filtered} nodes and {len(all_path_filtered)} pathways after filtration.')

            # filter invalid pathways
            filtered_pathways = reactions_filtration.filterPathways(tree_filtered)
            all_pathways_w_reactions = filtered_pathways

        # Check if we have at least 1 node and 1 pathway
        node_count = countNodes(tree_exp)
        all_path = searchPathways(tree_exp)

        if node_count < 1 or len(all_path) < 1:
            print(f"Warning: Insufficient data. The tree contains {node_count} nodes and {len(all_path)} pathways.")
            print(f"Saving raw results_dict data instead of pathways...")
            # Print the number of entries in results_dict
            print(f"Results dictionary contains {len(results_dict)} entries.")
            # Print a sample of the keys
            keys_sample = list(results_dict.keys())[:5] if len(results_dict) > 5 else list(results_dict.keys())
            print(f"Sample keys: {keys_sample}")
            # Use the raw results_dict instead of pathways when there are insufficient pathways
            all_pathways_w_reactions = results_dict

        # Format the pathways data properly
        # First, ensure we have a valid JSON structure
        formatted_pathways = {}

        # Check if all_pathways_w_reactions is already a dictionary
        if isinstance(all_pathways_w_reactions, dict):
            formatted_pathways = all_pathways_w_reactions
        else:
            # If it's not a dictionary, try to parse it
            try:
                # If it's a string that looks like a dictionary, try to parse it
                if isinstance(all_pathways_w_reactions, str) and all_pathways_w_reactions.strip().startswith('{'):
                    formatted_pathways = json.loads(all_pathways_w_reactions)
                else:
                    # Otherwise, just wrap it in a dictionary
                    formatted_pathways = {"pathways": all_pathways_w_reactions}
            except:
                # If parsing fails, just wrap it in a dictionary
                formatted_pathways = {"pathways": all_pathways_w_reactions}

        # Add metadata
        formatted_pathways["material"] = material
        formatted_pathways["retrieval_mode"] = retrieval_mode
        formatted_pathways["num_results"] = num_results
        formatted_pathways["alignment"] = alignment
        formatted_pathways["expansion"] = expansion
        formatted_pathways["filtration"] = filtration

        # Save the pathways to a JSON file
        print(f"Saving reaction pathways to {output_file}...")
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(formatted_pathways, f, indent=2, ensure_ascii=False)
            print(f"Successfully saved reaction pathways to {output_file}")
            # Print file size to confirm it was written
            file_size = os.path.getsize(output_file)
            print(f"File size: {file_size} bytes")
        except Exception as e:
            print(f"Error saving JSON file: {str(e)}")
            traceback.print_exc()

        # Return the pathways data
        return formatted_pathways

    except Exception as e:
        print(f"Error in json_main function: {str(e)}")
        traceback.print_exc()

        # Try to save error information to the output file
        try:
            error_data = {
                "error": str(e),
                "material": material,
                "retrieval_mode": retrieval_mode,
                "num_results": num_results,
                "alignment": alignment,
                "expansion": expansion,
                "filtration": filtration
            }

            print(f"Saving error information to {output_file}...")
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(error_data, f, indent=2, ensure_ascii=False)
            print(f"Successfully saved error information to {output_file}")
        except Exception as save_error:
            print(f"Failed to save error information: {str(save_error)}")

        return {"error": str(e)}


if __name__ == '__main__':
    try:
        args = parse_arguments()
        material = args.material
        num_results = args.num_results
        alignment = args.alignment == "True"
        expansion = args.expansion == "True"
        filtration = args.filtration == "True"
        retrieval_mode = args.retrieval_mode
        output_file = args.output

        print(
            f"Running with parameters: material={material}, num_results={num_results}, alignment={alignment}, expansion={expansion}, filtration={filtration}, retrieval_mode={retrieval_mode}, output={output_file}")

        result = main(
            material,
            num_results,
            alignment,
            expansion,
            filtration,
            retrieval_mode,
            output_file
        )
        print("Program completed successfully!")
    except Exception as e:
        print(f"Error in __main__: {str(e)}")
        traceback.print_exc()
