
"""
Base class for data generation: JSON reading, XML->JSON conversion.
"""
from pathlib import Path
import json
import xmltodict
from typing import Any

class DataGenerate:
    """Base class for data read/write operations."""
    
    def __init__(self, input_path: str | Path, output_path: str | Path = ""):
        self.input_path = Path(input_path) if input_path else None
        self.output_path = Path(output_path) if output_path else None
    
    def read_json(self, encoding: str = "utf-8") -> dict[str, Any]:
        """Reads a JSON file and returns a dictionary."""
        with open(self.input_path, 'r', encoding=encoding) as f:
            return json.load(f)
    
    def save_json(self, file_name: str | Path, data_dict: dict[str, Any]) -> None:
        """Saves a dictionary to a JSON file."""
        output = self.output_path / file_name if self.output_path else Path(file_name)
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, 'w', encoding='utf-8') as f:
            json.dump(data_dict, f, ensure_ascii=False, indent=2)
    
    def from_xml_to_json(self, encoding: str = "ISO-8859-1") -> dict[str, Any]:
        """Converts an XML file to a JSON dictionary."""
        with open(self.input_path, encoding=encoding) as xml_file:
            return xmltodict.parse(xml_file.read())