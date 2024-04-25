import json
import os

class Settings:
    """
    A class that represents the general settings for the 3D dust temperature map project, 
    by reading the general_settings.json file.
    It sets the locations of the temperature map data, code, plots, and paper location.
    """

    _instance = None

    def __new__(cls):
        """
        Creates a new instance of the Settings class if it doesn't already exist.

        Returns:
            The instance of the Settings class.
        """
        if cls._instance is None:
            cls._instance = super(Settings, cls).__new__(cls)
            try:
                settings_path = os.path.join(os.path.dirname(__file__), 'general_settings.json') # this adds the right slash for the OS that one is using
                with open(settings_path, 'r') as f:
                    cls._instance.data = json.load(f)
            except Exception as e:
                print(f"Failed to load configuration: {e}")
                exit(1)
        return cls._instance

    @property ## The property decorator allows the functions to be called without paranthesis, like an attribute, for better usage
    def data_location(self):
        
        return self.data.get('DataLocation', 'default/path/to/data') # the second parameter is the default value if the key is not found
    @property
    def code_location(self):
      
        return self.data.get('CodeLocation', 'default/path/to/code')
    
    @property
    def plots_location(self):
        return self.data.get('PlotsLocation', 'default/path/to/plots')
    @property
    def paper_location(self):
        return self.data.get('PaperLocation', 'default/path/to/paper')
    @property
    def configurations_location(self):
        return self.data.get('ConfigurationsLocation', os.path.join(os.path.dirname(__file__), 'configurations'))

###
# Usage in any module
# from general_settings import Settings

# settings = Settings()
# data_location = settings.temperature_map_data_location
# code_location = settings.temperature_map_code_location
