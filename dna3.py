import tkinter as tk
from tkinter import messagebox, scrolledtext
import requests
import matplotlib.pyplot as plt


def generate_rotations(sequence):
    return [sequence[i:] + sequence[:i] for i in range(len(sequence))]


def sort_rotations(rotations):
    return sorted(rotations)


def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))


def is_valid_dna(sequence):
    return all(base in 'ATCG' for base in sequence)


def sequence_statistics(sequence):
    length = len(sequence)
    gc_content = (sequence.count('G') + sequence.count('C')) / length * 100
    return f'Length: {length}, GC Content: {gc_content:.2f}%'


def find_pattern(sequence, pattern):
    occurrences = []
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i + len(pattern)] == pattern:
            occurrences.append(i)
    return occurrences

def extract_sequence_from_fasta(fasta_data):
    # Split the input data into lines
    lines = fasta_data.strip().split('\n')
    # The first line is usually the header with the name and description
    header = lines[0] if lines else ''
    # Join the lines that contain the sequence (which are not headers)
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    # Return both the header and the sequence
    return header, sequence



class MyApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DNA Sequence Visualizer")
        self.geometry("600x400")

        # Input text field
        self.input_text = tk.Entry(self, width=50)
        self.input_text.pack()

        # Process button
        self.process_button = tk.Button(self, text="GENERATE ROTATIONS", command=self.process_sequence)
        self.process_button.pack()

        # Output text area
        self.output_text = scrolledtext.ScrolledText(self, width=50, height=10)
        self.output_text.pack()

        # Reverse complement label
        self.reverse_comp_label = tk.Label(self, text="")
        self.reverse_comp_label.pack()

        # Sequence statistics label
        self.stats_label = tk.Label(self, text="")
        self.stats_label.pack()

        # Pattern search entry and button
        self.pattern_entry = tk.Entry(self, width=20)
        self.pattern_entry.pack()
        self.search_button = tk.Button(self, text="Search", command=self.search_pattern)
        self.search_button.pack()

        # Database selection and fetch button
        self.database_combobox = tk.StringVar(self)
        self.database_combobox.set("NCBI")
        tk.OptionMenu(self, self.database_combobox, "NCBI", "UniProt").pack()
        self.accession_entry = tk.Entry(self, width=20)
        self.accession_entry.pack()
        self.fetch_button = tk.Button(self, text="ANALYSE",command=lambda: self.fetch_sequence_from_database(self.accession_entry.get().strip(),self.database_combobox.get()))

        self.fetch_button.pack()

        # Visualize button
        self.visualize_button = tk.Button(self, text="VISUALIZE", command=self.visualize_sequence)
        self.visualize_button.pack()

    def process_sequence(self):
        sequence = self.input_text.get().strip().upper()
        if not is_valid_dna(sequence):
            messagebox.showwarning("Invalid Sequence", "The sequence contains invalid characters for DNA.")
            return
        rotations = generate_rotations(sequence)
        sorted_rotations = sort_rotations(rotations)
        self.display_rotations(sorted_rotations)
        rev_comp = reverse_complement(sequence)
        self.reverse_comp_label.config(text=f"Reverse Complement: {rev_comp}")
        stats = sequence_statistics(sequence)
        self.stats_label.config(text=stats)

    def display_rotations(self, rotations):
        self.output_text.delete(1.0, tk.END)
        for rotation in rotations:
            self.output_text.insert(tk.END, rotation + '\n')

    def search_pattern(self):
        sequence = self.input_text.get().strip().upper()
        pattern = self.pattern_entry.get().strip().upper()
        if not pattern:
            messagebox.showwarning("No Pattern", "Please enter a pattern to search for.")
            return
        occurrences = find_pattern(sequence, pattern)
        if occurrences:
            message = f'Pattern "{pattern}" found at positions: ' + ', '.join(str(pos) for pos in occurrences)
        else:
            message = f'Pattern "{pattern}" not found in the sequence.'
        messagebox.showinfo("Pattern Search Results", message)

    def fetch_sequence_from_database(self, accession_number, database_type):
        url = ''
        if database_type == "NCBI":
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession_number}&rettype=fasta&retmode=text"
        elif database_type == "UniProt":
            url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"

        try:
            response = requests.get(url)
            response.raise_for_status()  # Raises an HTTPError if the response was an HTTP error
            fasta_data = response.text
            header, sequence = extract_sequence_from_fasta(fasta_data)
            # You can now use 'header' to display the name
            messagebox.showinfo("Sequence Name", header)  # Display the name in a message box
            return sequence
        except requests.HTTPError as http_err:
            messagebox.showwarning("HTTP Error", f"HTTP error occurred: {http_err}")
        except requests.RequestException as err:
            messagebox.showwarning("Error Fetching Sequence", f"An error occurred: {err}")
        return "", ""

    def visualize_sequence(self):
        patient_sequence = self.input_text.get().strip().upper()

        # Fetch the reference sequence using the accession number and database selection
        accession_number = self.accession_entry.get().strip()
        database_type = self.database_combobox.get()
        reference_sequence = self.fetch_sequence_from_database(accession_number, database_type)

        # Calculate nucleotide composition for both sequences
        patient_counts = {'A': patient_sequence.count('A'), 'T': patient_sequence.count('T'),
                          'G': patient_sequence.count('G'), 'C': patient_sequence.count('C')}
        reference_counts = {'A': reference_sequence.count('A'), 'T': reference_sequence.count('T'),
                            'G': reference_sequence.count('G'), 'C': reference_sequence.count('C')}


        nucleotides = ['A', 'T', 'G', 'C']
        colors = ['red', 'blue', 'green', 'yellow']

        # Normalize counts to proportions for pie charts
        total_patient = sum(patient_counts.values())
        patient_proportions = [count / total_patient for count in patient_counts.values()]

        total_reference = sum(reference_counts.values())
        reference_proportions = [count / total_reference for count in reference_counts.values()]

        explode_values = (0.1, 0, 0, 0)  # Only the first slice (A) is popped out

        # Create a figure with two subplots for pie charts
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

        # Plot pie chart for patient sequence


        ax1.pie(patient_proportions, labels=nucleotides, colors=colors, autopct='%1.1f%%', startangle=90, explode=explode_values, shadow=True)
        ax1.set_title('Provided DNA Composition', fontsize=14)

        # Plot pie chart for reference sequence
        ax2.pie(reference_proportions, labels=nucleotides, colors=colors, autopct='%1.1f%%', startangle=90,
                explode=explode_values, shadow=True)
        ax2.set_title("A Diseased patient's DNA composition", fontsize=14)

        # Add a legend
        fig.legend(nucleotides, loc='lower center', ncol=4, fontsize=12)
       # Adjust layout to make room for the legend
        plt.subplots_adjust(bottom=0.2)

        plt.tight_layout()
        plt.show()




if __name__ == "__main__":
    app = MyApp()
    app.mainloop()
