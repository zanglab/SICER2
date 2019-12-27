
class Sicer():
    def __init__(
        self,
        treatment_file,
        control_file,
        species,
        redundancy_threshold,
        window_size,
        fragment_size,
        effective_genome_fraction,
        false_discovery_rate,
        output_directory,
        gap_size,
        e_value,
        cpu,
        significant_reads,
        verbose,
        df):

            self.treatment_file = treatment_file
            self.control_file = control_file
            self.species = species
            self.redundancy_threshold = redundancy_threshold
            self.window_size = window_size
            self.fragment_size = fragment_size
            self.effective_genome_fraction = effective_genome_fraction
            self.false_discovery_rate = false_discovery_rate
            self.output_directory = output_directory
            self.gap_size = gap_size
            self.e_value = e_value
            self.cpu = cpu
            self.significant_reads = significant_reads
            self.verbose = verbose
            self.df = df


    def run():
        

