from Bio import AlignIO
from Bio import SeqIO
import random
import pickle



def parameter_calculation(year_n_fasta, year_n_plus_one_fasta):
    """
    Parameter_Calculation does multiple sequence alignment then

    :param year_n_fasta: file name containing cRNA sequence of gene(s)
    :param year_n_plus_one_fasta: file name containing cRNA sequence of same gene(s) of the next year
    :return:
    """
    codon_variation_dictionary = {}
    codon_transition_dictionary = {}
    total_nucleotides = 0
    with open(year_n_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            total_nucleotides = total_nucleotides + len(record.seq)
            with open("alignment.fasta", "w") as file_handle:
                SeqIO.write(record, file_handle, "fasta")

                for other_records in SeqIO.parse(year_n_plus_one_fasta, "fasta"):
                    total_nucleotides = total_nucleotides + len(record.seq)
                    if len(record.seq) == len(other_records.seq):
                        SeqIO.write(other_records, file_handle, "fasta")
            file_handle.close()
            """
            Multiple Sequence Alignment
            """
            for multiple_alignment in AlignIO.parse("alignment.fasta", "fasta"):
                i = 0
                base_seq = []
                for seqreq in multiple_alignment:
                    if i == 0:
                        base_seq = seqreq
                        i += 1
                    else:
                        i += 1
                        for j in range(0, len(seqreq.seq) - 3, 3):
                            base_codon = str(base_seq.seq[j:j + 3])
                            compared_codon = str(seqreq.seq[j:j + 3])

                            base_codon = base_codon.replace("R",random.choice(["A","U"]))
                            base_codon = base_codon.replace("T","U")
                            base_codon = base_codon.replace("Y",random.choice(["C","U"]))
                            compared_codon = compared_codon.replace("R",random.choice(["A","U"]))
                            compared_codon = compared_codon.replace("T","U")
                            compared_codon = compared_codon.replace("Y",random.choice(["C","U"]))

                            if codon_variation_dictionary.get(base_codon,0)==0:
                                codon_variation_dictionary[base_codon] ={}
                            else:
                                codon_variation_dictionary[base_codon][compared_codon] = codon_variation_dictionary.get(base_codon, {}).get(compared_codon, 0) + 1
    """
    Transition Matrix conversion
    """
    for row_codon in codon_variation_dictionary.keys():
        total_transitions = 0
        for col_codon in codon_variation_dictionary[row_codon].keys():
            total_transitions += codon_variation_dictionary.get(row_codon, 0).get(col_codon, 0)
        for col_codon in codon_variation_dictionary[row_codon].keys():
            value = codon_variation_dictionary.get(row_codon, 0).get(col_codon, 0)
            codon_transition_dictionary[row_codon] = codon_transition_dictionary.get(row_codon,{})
            codon_transition_dictionary[row_codon][col_codon] = value / total_transitions
    print(codon_transition_dictionary)

    with open("transition.pkl","wb") as pkl:
        pickle.dump(codon_transition_dictionary,pkl)


def random_sequence_generation(seed_sequence_fasta, num_seq):
    """

    :param seed_sequence_fasta:
    :param num_seq:
    :return:
    """

    with open("transition.pkl","rb") as pkl:
        transition_dictionary = pickle.load(pkl)

    with open("generated_sequences.fasta", "w") as file_handle:
        for seq in SeqIO.parse(seed_sequence_fasta,"fasta"):
            for j in range(num_seq):
                gen_seq = ""
                for i in range(0,len(seq.seq)-3,3):
                    codon = str(seq.seq[i:i+3])

                    codon = codon.replace("R", random.choice(["A", "U"]))
                    codon = codon.replace("T", "U")
                    codon = codon.replace("Y", random.choice(["C", "U"]))
                    markov_states = []
                    markov_weights = []
                    for col_codon in transition_dictionary[codon].keys():
                        markov_states.append(col_codon)
                        markov_weights.append(transition_dictionary[codon][col_codon])
                    print(random.choices(markov_states,markov_weights)[0])
                    gen_seq = gen_seq+random.choices(markov_states,markov_weights)[0]
                with open("generated_sequences.fasta","a") as fasta_writer:
                    fasta_writer.write(">Generated Sequence:"+str(j)+"\n")
                    fasta_writer.write(gen_seq+"\n")









def validation():
    return "Stub"

from Bio import AlignIO
from Bio import SeqIO
import random
import pickle



def parameter_calculation(year_n_fasta, year_n_plus_one_fasta):
    """
    Parameter_Calculation does multiple sequence alignment then

    :param year_n_fasta: file name containing cRNA sequence of gene(s)
    :param year_n_plus_one_fasta: file name containing cRNA sequence of same gene(s) of the next year
    :return:
    """
    codon_variation_dictionary = {}
    codon_transition_dictionary = {}
    total_nucleotides = 0
    with open(year_n_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            total_nucleotides = total_nucleotides + len(record.seq)
            with open("alignment.fasta", "w") as file_handle:
                SeqIO.write(record, file_handle, "fasta")

                for other_records in SeqIO.parse(year_n_plus_one_fasta, "fasta"):
                    total_nucleotides = total_nucleotides + len(record.seq)
                    if len(record.seq) == len(other_records.seq):
                        SeqIO.write(other_records, file_handle, "fasta")
            file_handle.close()
            """
            Multiple Sequence Alignment
            """
            for multiple_alignment in AlignIO.parse("alignment.fasta", "fasta"):
                i = 0
                base_seq = []
                for seqreq in multiple_alignment:
                    if i == 0:
                        base_seq = seqreq
                        i += 1
                    else:
                        i += 1
                        for j in range(0, len(seqreq.seq) - 3, 3):
                            base_codon = str(base_seq.seq[j:j + 3])
                            compared_codon = str(seqreq.seq[j:j + 3])

                            base_codon = base_codon.replace("R",random.choice(["A","U"]))
                            base_codon = base_codon.replace("T","U")
                            base_codon = base_codon.replace("Y",random.choice(["C","U"]))
                            compared_codon = compared_codon.replace("R",random.choice(["A","U"]))
                            compared_codon = compared_codon.replace("T","U")
                            compared_codon = compared_codon.replace("Y",random.choice(["C","U"]))

                            if codon_variation_dictionary.get(base_codon,0)==0:
                                codon_variation_dictionary[base_codon] ={}
                            else:
                                codon_variation_dictionary[base_codon][compared_codon] = codon_variation_dictionary.get(base_codon, {}).get(compared_codon, 0) + 1
    """
    Transition Matrix conversion
    """
    for row_codon in codon_variation_dictionary.keys():
        total_transitions = 0
        for col_codon in codon_variation_dictionary[row_codon].keys():
            total_transitions += codon_variation_dictionary.get(row_codon, 0).get(col_codon, 0)
        for col_codon in codon_variation_dictionary[row_codon].keys():
            value = codon_variation_dictionary.get(row_codon, 0).get(col_codon, 0)
            codon_transition_dictionary[row_codon] = codon_transition_dictionary.get(row_codon,{})
            codon_transition_dictionary[row_codon][col_codon] = value / total_transitions
    print(codon_transition_dictionary)

    with open("transition.pkl","wb") as pkl:
        pickle.dump(codon_transition_dictionary,pkl)


def random_sequence_generation(seed_sequence_fasta, num_seq):
    """

    :param seed_sequence_fasta:
    :param num_seq:
    :return:
    """

    with open("transition.pkl","rb") as pkl:
        transition_dictionary = pickle.load(pkl)

    with open("generated_sequences.fasta", "w") as file_handle:
        for seq in SeqIO.parse(seed_sequence_fasta,"fasta"):
            for j in range(num_seq):
                gen_seq = ""
                for i in range(0,len(seq.seq)-3,3):
                    codon = str(seq.seq[i:i+3])

                    codon = codon.replace("R", random.choice(["A", "U"]))
                    codon = codon.replace("T", "U")
                    codon = codon.replace("Y", random.choice(["C", "U"]))
                    markov_states = []
                    markov_weights = []
                    for col_codon in transition_dictionary[codon].keys():
                        markov_states.append(col_codon)
                        markov_weights.append(transition_dictionary[codon][col_codon])
                    print(random.choices(markov_states,markov_weights)[0])
                    gen_seq = gen_seq+random.choices(markov_states,markov_weights)[0]
                with open("generated_sequences.fasta","a") as fasta_writer:
                    fasta_writer.write(">Generated Sequence:"+str(j)+"\n")
                    fasta_writer.write(gen_seq+"\n")









def validation():
    return "Stub"

