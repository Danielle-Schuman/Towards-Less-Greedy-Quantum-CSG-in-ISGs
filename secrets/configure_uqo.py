# Step 1: run this to create your uqo private key file
import os
from uqo.generate_certificates import generate_certificates

generate_certificates(os.path.dirname(os.path.realpath(__file__)))

# Step 2: Create a file named "config.json" with the following content:
# {
#      "method": "token",
#      "endpoint": "uq.mobile.ifi.lmu.de:30000",
#      "credentials": "YOUR_TOKEN",
#      "private_key_file": "secrets/private_keys/client.key_secret"
# }
#
# where YOUR_TOKEN is your uqo-token