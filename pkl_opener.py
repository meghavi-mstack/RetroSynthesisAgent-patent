import pickle
import json

# Load the pickle file
with open('tree_pi/Azulene_w_exp_alg.pkl', 'rb') as file:
    my_obj = pickle.load(file)

# Function to convert any object to JSON-serializable dict
def safe_serialize(obj):
    if isinstance(obj, (str, int, float, bool, type(None))):
        return obj
    elif isinstance(obj, dict):
        return {str(k): safe_serialize(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple, set)):
        return [safe_serialize(i) for i in obj]
    elif hasattr(obj, '__dict__'):
        return safe_serialize(obj.__dict__)
    else:
        return str(obj)  # fallback to string representation

# Build a dict from dir() and values
data = {}
for attr in dir(my_obj):
    if attr.startswith('__') and attr.endswith('__'):
        continue  # Skip dunder methods
    try:
        value = getattr(my_obj, attr)
        data[attr] = safe_serialize(value)
    except Exception as e:
        data[attr] = f"Could not access attribute: {str(e)}"

# Write to JSON file
with open('tree_object_dump.json', 'w') as f:
    json.dump(data, f, indent=2)

print("✅ Tree object data dumped into tree_object_dump.json")
