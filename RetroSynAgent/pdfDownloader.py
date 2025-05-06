import os
import re
import requests
from loguru import logger
import json
import random
import glob
import difflib
from concurrent.futures import ThreadPoolExecutor, as_completed
from scholarly import scholarly
from dotenv import load_dotenv
import os

class PDFDownloader:
    def __init__(self, material, pdf_folder_name, num_results = 5, n_thread=3):
        self.pdf_folder_name = pdf_folder_name
        self.no_download_link_json_name = 'no_download_link_titles.json'
        self.query = material + ' AND synthesis'
        self.num_results = num_results  # Store num_results for later use
        self.title_list = self.get_scholar_titles(self.query, num_results)
        # self.title_list = self.get_scholar_dois(self.query, num_results)
        self.n_thread = n_thread
        self.url = 'https://www.sci-hub.se/'
        self.no_download_link_titles = self.read_data_from_json(self.no_download_link_json_name) if os.path.exists(
            self.no_download_link_json_name) else []

        load_dotenv()
        headers_dict = os.getenv("HEADERS")
        self.headers = json.loads(headers_dict)
        cookies_dict = os.getenv("COOKIES")
        self.cookies = json.loads(cookies_dict)

    # version2
    def get_scholar_titles(self, query, num_results, citations=0):
        import time
        import threading

        title_list = []
        search_completed = False

        def search_with_timeout():
            nonlocal title_list, search_completed
            try:
                search_query = scholarly.search_pubs(query)
                count = 0
                for i in range(10000):
                    try:
                        pub = next(search_query)
                        if pub['num_citations'] > citations:
                            title = pub['bib']['title']
                            num_citations = int(pub['num_citations'])
                            title_list.append(title)
                            count += 1
                            if count == num_results:
                                break
                    except StopIteration:
                        break
                    except KeyError:
                        continue
                    except Exception as e:
                        print(f"Error processing publication: {str(e)}")
                        continue
                search_completed = True
            except Exception as e:
                print(f"Error in Google Scholar search: {str(e)}")
                search_completed = True

        # Start the search in a separate thread
        search_thread = threading.Thread(target=search_with_timeout)
        search_thread.daemon = True
        search_thread.start()

        # Wait for up to 60 seconds
        timeout = 60  # seconds
        start_time = time.time()
        while not search_completed and time.time() - start_time < timeout:
            time.sleep(1)

        if not search_completed:
            print(f"Google Scholar search timed out after {timeout} seconds")
            print("Using empty title list and proceeding with alternative sources")

        return title_list

    def read_data_from_json(self, filename):
        with open(filename, 'r', encoding='utf-8') as file:
            data = json.load(file)
        return data

    def save_data_as_json(self, filename, data):
        with open(filename, 'w') as json_file:
            json.dump(data, json_file, indent=4, ensure_ascii=False)
        # print(f'{filename} saved!')

    # def get_pdf_files(self):
    #     pdf_files = glob.glob(os.path.join(self.pdf_folder_name, '*.pdf'))
    #     return [os.path.basename(file) for file in pdf_files]

    def get_pdf_files(self):
        # get pdf path list (pdf/pdf_sub/title.pdf)
        # pdf_files = []
        # for root, dirs, files in os.walk(self.pdf_folder_name):
        #     for file in files:
        #         if file.endswith(".pdf"):
        #             pdf_files.append(os.path.join(root, file))
        # return pdf_files
        # get pdf name list (title.pdf)
        pdf_files = []
        for root, dirs, files in os.walk(self.pdf_folder_name):
            for file in files:
                if file.endswith(".pdf"):
                    pdf_files.append(file)  # 只获取文件名
        return pdf_files

    def check_pdf_existence(self, target_pdf_name, pdf_name_list):
        similarity_threshold = 0.9
        for pdf_name in pdf_name_list:
            similarity = difflib.SequenceMatcher(None, target_pdf_name, pdf_name).ratio()
            if similarity > similarity_threshold:
                return True
        return False

    def title_href(self, title):
        data = {"request": str(title)}
        try:
            res = requests.post(self.url, headers=self.headers, data=data, cookies=self.cookies,verify=False)
            _href = re.findall("""<button onclick = "location.href='(.*?)'">&darr;""", res.text)
            if _href:
                _href = _href[0]
                if 'sci-hub' not in _href:
                    download_href = "https://sci-hub.se" + _href
                else:
                    download_href = "https:" + _href

                return download_href

            elif '<p id = "smile">:(</p>' in res.text:
                # logger.error(f"No download link for {title} in sci-hub")
                self.no_download_link_titles.append(title)
                self.save_data_as_json(self.no_download_link_json_name, self.no_download_link_titles)
                return None
            else:
                # logger.error(f"Failed to get download link for title: {title}, status_code: {res.status_code}")
                return None
        except Exception as e:
            # logger.error(e)
            return None

    def get_download_pdf(self, href, title):
        try:
            res = requests.get(href, headers=self.headers, cookies=self.cookies,verify=False)
            if res.status_code == 200:
                file_name = f"{str(title).replace(':', '').replace('/', '').replace('*', '').replace('|', '').replace('?', '')}.pdf"
                file_path = os.path.join(os.getcwd(), self.pdf_folder_name, file_name)
                if res.headers.get('Content-Type') == 'application/pdf':
                    with open(file_path, 'wb') as f:
                        f.write(res.content)
                        logger.info(f"{file_name} successfully saved!")
                else:
                    logger.error(f"{file_name} is invalid pdf file")
            else:
                logger.error(f"Failed to download pdf for {title}, status code: {res.status_code}")
                if res.status_code == 404:
                    self.no_download_link_titles.append(title)
                    self.save_data_as_json(self.no_download_link_json_name, self.no_download_link_titles)
        except Exception as e:
            logger.error(e)

    def download_task(self, title):
        href = self.title_href(title)
        if href:
            self.get_download_pdf(href, title)

    def download_pdfs(self, titles):
        with ThreadPoolExecutor(max_workers=self.n_thread) as executor:
            futures = [executor.submit(self.download_task, title) for title in titles]
            for future in as_completed(futures):
                future.result()

    def filter_titles(self, titles):
        pdf_file_list = self.get_pdf_files()
        pdf_name_list = [pdf.split('.pdf')[0] for pdf in pdf_file_list]

        titles_filtered1 = [title for title in titles if not self.check_pdf_existence(title, pdf_name_list)]
        titles_filtered2 = [title for title in titles_filtered1 if
                            not self.check_pdf_existence(title, self.no_download_link_titles)]
        return titles_filtered2

    def get_fallback_titles(self, query, num_results=10):
        """
        Fallback method to generate titles when Google Scholar fails
        """
        # Extract the main substance name from the query (remove 'AND synthesis')
        substance = query.replace(' AND synthesis', '')

        # Generate some generic titles related to the substance
        generic_titles = [
            f"Synthesis and Characterization of {substance}",
            f"Novel Synthetic Routes to {substance}",
            f"Efficient Synthesis of {substance} and Its Derivatives",
            f"A Review of Synthetic Approaches to {substance}",
            f"Green Synthesis of {substance}",
            f"One-Pot Synthesis of {substance}",
            f"Catalytic Synthesis of {substance}",
            f"Improved Method for the Preparation of {substance}",
            f"Recent Advances in the Synthesis of {substance}",
            f"Sustainable Synthesis of {substance}"
        ]

        # Return only the requested number of titles
        return generic_titles[:num_results]

    def main(self):
        print(f"PDFDownloader: Creating directory {self.pdf_folder_name}")
        os.makedirs(self.pdf_folder_name, exist_ok=True)
        print(f"PDFDownloader: Query = {self.query}")

        # Check if we got any titles from Google Scholar
        if not self.title_list:
            print("PDFDownloader: No titles found from Google Scholar. Using fallback titles.")
            self.title_list = self.get_fallback_titles(self.query, self.num_results)
            print(f"PDFDownloader: Generated {len(self.title_list)} fallback titles")
        else:
            print(f"PDFDownloader: Found {len(self.title_list)} titles from Google Scholar")

        joined_titles = '\n'.join(self.title_list)
        print(f"PDFDownloader: Title list:\n{joined_titles}")

        print("PDFDownloader: Filtering titles...")
        titles_filtered = self.filter_titles(self.title_list)

        print(f'PDFDownloader: Total number of titles: {len(self.title_list)}, '
              f'{len(self.get_pdf_files())} have been downloaded, '
              f'{len(self.no_download_link_titles)} do not have download links, '
              f'{len(titles_filtered)} are planned to be downloaded.')

        # random.shuffle(titles_filtered)
        print("PDFDownloader: Starting download of PDFs...")
        self.download_pdfs(titles_filtered)

        download_pdf_filename_list = self.get_pdf_files()
        print(f'PDFDownloader: Finished downloading, '
              f'{len(download_pdf_filename_list)} have been downloaded, '
              f'{len(self.no_download_link_titles)} do not have download links, '
              f'{len(self.title_list)-len(self.get_pdf_files())-len(self.no_download_link_titles)} failed to be downloaded.')

        # Return full paths to the downloaded PDFs
        full_paths = []
        for pdf_name in download_pdf_filename_list:
            full_path = os.path.join(os.path.abspath(self.pdf_folder_name), pdf_name)
            full_paths.append(full_path)
        return full_paths
