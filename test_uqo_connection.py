from uqo.client.config import Config
from uqo.examples.examples import *

config = Config(configpath="secrets_folder/config.json")
connection = config.create_connection()

ping(config)
print("Available solvers: ", connection.get_available_dwave_solvers())
print("Remaining Quota: ")
connection.show_quota()
