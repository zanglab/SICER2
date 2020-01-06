import os
import json

genomedata_dir = os.path.dirname(os.path.dirname(__file__)) + "/genomedata"

def get_available_species():
    available_species = []
    for file in os.listdir(genomedata_dir):
        if file.endswith(".json"):
            available_species.append(os.path.splitext(file)[0])
    return available_species


class GenomeData:
    # Class used to handle chromsome meta-data

    def __init__(self, species):

        if species in get_available_species():
            self.species = species
            file = genomedata_dir + "/" + species + ".json"
        else:
            raise ValueError("\"" + species + "\" is not an available species.")

        with open(file) as fp:
            data = json.load(fp)
            self.chrom = data["chrom"]
            self.chrom_length = data["chrom_length"]

