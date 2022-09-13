from pathlib import Path
import json

file_paths = sorted(Path('.').glob('hello_world[1-2].txt'))

def add_data(paths):
    data = []
    for path in paths:
        data.append(
            {
                "HelloWorld.hello": f"{path}"
            }
        )
    return data


if __name__=="__main__":
    output = Path('hello_world_batch.json')
    with output.open('w') as f:
        json.dump(add_data(file_paths),f,indent=4)  

