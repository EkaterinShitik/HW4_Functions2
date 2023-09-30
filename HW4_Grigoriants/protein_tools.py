import dictionaries


def three_one_letter_code(sequences: str) -> list:
    """
    Reverse the protein sequences from one-letter to three-letter format and vice-versa

    Case 1: get three-letter sequence\n
    Use one-letter amino-acids sequences of any letter case

    Case 2: get one-letter sequence\n
    Use three-letter amino-acid separated by "-" sequences.
    Please note that sequences without "-" are parsed as one-letter code sequences\n
    Example: for sequence "Ala" function will return "Ala-leu-ala"

    Arguments:
    - sequences (tuple[str] or list[str]): protein sequences to convert\n
    Example: ["WAG", "MkqRe", "msrlk", "Met-Ala-Gly", "Met-arg-asn-Trp-Ala-Gly", "arg-asn-trp"]

    Return:
    - list: one-letter/three-letter protein sequences\n
    Example: ["Met-Ala-Gly", "Met-arg-asn-Trp-Ala-Gly", "arg-asn-trp", "WAG", "MkqRe", "rlk"]
    """
    inversed_sequences = []
    for sequence in sequences:
        inversed_sequence = ""
        if "-" not in sequence:
            for letter in sequence:
                inversed_sequence += dictionaries.amino_acids[letter] + "-"
            inversed_sequence = inversed_sequence[:-1]
            inversed_sequences.append(inversed_sequence)
        else:
            aa_splitted = sequence.split("-")
            for aa in aa_splitted:
                inversed_sequence += list(dictionaries.amino_acids.keys())[
                    list(dictionaries.amino_acids.values()).index(aa)
                ]
            inversed_sequences.append(inversed_sequence)
    return inversed_sequences


def define_molecular_weight(sequences: str) -> dict:
    """
    Define molecular weight of the protein sequences

    Use one-letter amino-acids sequences of any letter case
    The molecular weight is:
    - a sum of masses of each atom constituting a molecule
    - expressed in units called daltons (Da)
    - rounded to hundredths

    Arguments:
    - sequences (tuple[str] or list[str]): protein sequences to convert

    Return:
    - dictionary: protein sequences as keys and molecular masses as values\n
    Example: {"WAG": 332.39, "MkqRe": 690.88, "msrlk": 633.86}
    """
    sequences_weights = {}
    for sequence in sequences:
        sequence_weight = 0
        for letter in sequence:
            sequence_weight += dictionaries.amino_acid_weights[letter.upper()]
        sequence_weight -= (len(sequence) - 1) * 18  # deduct water from peptide bond
        sequences_weights[sequence] = round(sequence_weight, 2)
    return sequences_weights


def search_for_motifs(
    sequences: (tuple[str] or list[str]), motif: str, overlapping: bool
) -> dict:
    """
    Search for motifs - conserved amino acids residues in protein sequence

    Search for one motif at a time\n
    Search is letter case sensitive\n
    Use one-letter aminoacids code for desired sequences and motifs\n
    Positions of AA in sequences are counted from 0\n
    By default, overlapping matches are counted

    Arguments:
    - sequences (tuple[str] or list[str]): sequences to check for given motif within\n
        Example: sequences = ["AMGAGW", "GAWSGRAGA"]
    - motif (str]: desired motif to check presense in every given sequence\n
        Example: motif = "GA"
    - overlapping (bool): count (True) or skip (False) overlapping matches. (Optional)\n
        Example: overlapping = False
    Return:
    - dictionary: sequences (str] as keys , starting positions for presented motif (list) as values\n
        Example: {"AMGAGW": [2], "GAWSGRAGA": [0, 7]}
    """
    new_line = "\n"
    all_positions = {}
    for sequence in sequences:
        start = 0
        positions = []
        print(f"Sequence: {sequence}")
        print(f"Motif: {motif}")
        if motif in sequence:
            while True:
                start = sequence.find(motif, start)
                if start == -1:
                    break
                positions.append(start)
                if overlapping:
                    start += 1
                else:
                    start += len(motif)
            print_pos = ", ".join(str(x) for x in positions)
            print_pos = f"{print_pos}{new_line}"
            print(
                f"Motif is present in protein sequence starting at positions: {print_pos}"
            )
        else:
            print(f"Motif is not present in protein sequence{new_line}")
        all_positions[sequence] = positions
    return all_positions


def search_for_alt_frames(sequences: str, alt_start_aa: str) -> dict:
    """
    Search for alternative frames in a protein sequences

    Search is not letter case sensitive\n
    Without an alt_start_aa argument search for frames that start with methionine ("M")
    To search frames with alternative start codon add alt_start_aa argument\n
    In alt_start_aa argument use one-letter code

    The function ignores the last three amino acids in sequences

    Arguments:
    - sequences (tuple[str] or list[str]): sequences to check
    - alt_start_aa (str]: the name of an amino acid that is encoded by alternative start AA (Optional)\n
    Example: alt_start_aa = "I"

    Return:
    - dictionary: the original sequences and a collection of alternative frames
    """
    alternative_frames = {}
    num_position = 0
    for sequence in sequences:
        alternative_frames[sequence] = []
        for amino_acid in sequence[1:-3]:
            alt_frame = ""
            num_position += 1
            if amino_acid == alt_start_aa or amino_acid == alt_start_aa.swapcase():
                alt_frame += sequence[num_position:]
                alternative_frames[sequence].append(alt_frame)
        num_position = 0
    return alternative_frames


def convert_to_nucl_acids(sequences: list, nucl_acids: str) -> dict:
    """
    Convert protein sequences to RNA or DNA sequences.

    Use the most frequent codons in human. The source - https://www.genscript.com/tools/codon-frequency-table\n
    All nucleic acids (DNA and RNA) are showed in 5"-3" direction

    Arguments:
    - sequences (tuple[str] or list[str]): sequences to convert
    - nucl_acids (str]: the nucleic acid that is prefered\n
    Example: nucl_acids = "RNA" - convert to RNA\n
                   nucl_acids = "DNA" - convert to DNA\n
                   nucl_acids = "both" - convert to RNA and DNA
    Return:
    - dictionary: nucleic acids (str) as keys, collection of sequences (list) as values
    """
    rule_of_translation = sequences[0].maketrans(dictionaries.translation_rule)
    rule_of_transcription = sequences[0].maketrans("AaUuCcGg", "TtAaGgCc")
    nucl_acid_seqs = {"RNA": [], "DNA": []}
    for sequence in sequences:
        rna_seq = sequence.translate(rule_of_translation)
        dna_seq = rna_seq.translate(rule_of_transcription)
        if nucl_acids == "RNA":
            nucl_acid_seqs["RNA"].append(rna_seq)
            if sequence == sequences[-1]:
                del nucl_acid_seqs["DNA"]
        if nucl_acids == "DNA":
            nucl_acid_seqs["DNA"].append(dna_seq)
            if sequence == sequences[-1]:
                del nucl_acid_seqs["RNA"]
        if nucl_acids == "both":
            nucl_acid_seqs["RNA"].append(rna_seq)
            nucl_acid_seqs["DNA"].append(dna_seq)
    return nucl_acid_seqs


procedures_to_functions = {
    "search_for_motifs": search_for_motifs,
    "search_for_alt_frames": search_for_alt_frames,
    "convert_to_nucl_acids": convert_to_nucl_acids,
    "three_one_letter_code": three_one_letter_code,
    "define_molecular_weight": define_molecular_weight,
}


def check_and_parse_user_input(
    sequences: list[str] or tuple[str], **kwargs
) -> dict and str:
    """
    Check if user input can be correctly processed\n
    Parse sequences and arguments for desired procedure
    
    Arguments:
    - sequences (list[str] or tuple[str]): sequences to process
    - **kwargs - needed arguments for completion of desired procedure

    Return:
    - string: procedure name
    - dictionary: a collection of procedure arguments and their values
    """
    if len(sequences) == 0:
        raise ValueError("No sequences provided")
    procedure = kwargs["procedure"]
    if procedure not in procedures_to_functions.keys():
        raise ValueError("Wrong procedure")
    allowed_inputs = set(dictionaries.amino_acids.keys()).union(
        set(dictionaries.amino_acids.values())
    )
    allowed_inputs.add("-")
    if procedure != "three_one_letter_code":
        allowed_inputs -= set(dictionaries.amino_acids.values())
    for sequence in sequences:
        allowed_inputs_seq = allowed_inputs.copy()
        if procedure == "three_one_letter_code" and "-" in sequence:
            allowed_inputs_seq -= set(dictionaries.amino_acids.keys())
            if not all(
                aminoacids in allowed_inputs_seq for aminoacids in sequence.split("-")
            ):
                raise ValueError("Invalid sequence given")
        else:
            allowed_inputs_seq.remove("-")
            allowed_inputs_seq -= set(dictionaries.amino_acids.values())
            if not all(aminoacids in allowed_inputs_seq for aminoacids in sequence):
                raise ValueError("Invalid sequence given")
    procedure_arguments = {}
    if procedure == "search_for_motifs":
        if "motif" not in kwargs.keys():
            raise ValueError("Please provide desired motif")
        procedure_arguments["motif"] = kwargs["motif"]
        if "overlapping" not in kwargs.keys():
            procedure_arguments["overlapping"] = True
        else:
            procedure_arguments["overlapping"] = kwargs["overlapping"]
    elif procedure == "search_for_alt_frames":
        if "alt_start_aa" not in kwargs.keys():
            procedure_arguments["alt_start_aa"] = "M"
        else:
            if len(kwargs["alt_start_aa"]) > 1:
                raise ValueError("Invalid alternative start AA")
            procedure_arguments["alt_start_aa"] = kwargs["alt_start_aa"]
    elif procedure == "convert_to_nucl_acids":
        if "nucl_acids" not in kwargs.keys():
            raise ValueError("Please provide desired type of nucl_acids")
        if kwargs["nucl_acids"] not in {"DNA", "RNA", "both"}:
            raise ValueError("Invalid nucl_acids argument")
        procedure_arguments["nucl_acids"] = kwargs["nucl_acids"]
    procedure_arguments["sequences"] = sequences
    return procedure_arguments, procedure


def run_protein_tools(sequences: list[str] or tuple[str], **kwargs: str):
    """
    Main function to process protein sequence by one of the developed tools.\n
    Run one procedure at a time:
    - Search for conserved amino acids residues in protein sequence
    - Search for alternative frames in a protein sequences
    - Convert protein sequences to RNA or DNA sequences
    - Reverse the protein sequences from one-letter to three-letter format and vice-versa
    - Define molecular weight of the protein sequences

    All functions except *search_for_alt_frames* are letter case sensitive\n
    Provide protein sequence in one letter code.\n
    You can obtain one letter code from three letter code with *three_one_letter_code*\n
    If more information needed please see README or desired docstring

    Arguments:
    - sequences (list[str] or tuple[str]): sequences to process
    - procedure (str]: desired procedure:
        - "search_for_motifs"
        - "search_for_alt_frames"
        - "convert_to_nucl_acids"
        - "three_one_letter_code"
        - "define_molecular_weight"

    For "search_for_motif" procedure provide:
    - motif (str]: desired motif to check presense in every given sequence\n
            Example: motif = "GA"
    - overlapping (bool): count (True) or skip (False) overlapping matches. (Optional)\n
            Example: overlapping = False

    For "search_for_alt_frames" procedure provide:
    - alt_start_aa (str]: the name of an amino acid that is encoded by alternative start codon (Optional)\n
            Example: alt_start_aa = "I"

    For "convert_to_nucl_acids" procedure provide:
    - nucl_acids (str]: the nucleic acid to convert to\n
            Example: nucl_acids = "RNA"\n
                           nucl_acids = "DNA"\n
                           nucl_acids = "both"

    Return:
    - dict: Dictionary with processed sequences. Depends on desired tool\n
            Please see Readme or desired docstring
    """
    procedure_arguments, procedure = check_and_parse_user_input(sequences, **kwargs)
    return procedures_to_functions[procedure](**procedure_arguments)
