import os,io
from Bio import SeqIO
from tqdm import tqdm

def formatfasta(ana, seq, num=70)->str:
    """
    format sequence from string to fasta
    :param ana: analysis
    :param seq: sequence
    :param num: number of characters per line
    :return: fasta string
    """
    ana = ana.replace('\n','').replace('\r','')

    format_seq = ""
    for i, char in enumerate(seq):
        format_seq += char
        if (i + 1) % num == 0:
            format_seq += "\n"
    return ana +'\n' + format_seq + "\n"


def splitfasta(fasta_file, out_dir)->int:
    """
    split fasta file into multiple fasta files
    :param fasta_file: fasta file
    :param out_dir: output directory
    :return: sequence count
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    seq_num = 0
    with io.open(fasta_file) as f:
        fileline = ""
        for line in tqdm(f.readlines()):
            if line.startswith(">"):
                seq_num += 1
                if fileline:
                    with io.open(os.path.join(out_dir, "seq_%d.fasta" % seq_num), "w") as f:
                        f.write(fileline)
                    fileline = ""
            fileline+=line
    return seq_num


def splitgb(gb_file, out_dir)->int:
    """
    split genbank file into multiple genbank files
    :param gb_file: genbank file
    :param out_dir: output directory
    :return: sequence count
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    seq_num = 0
    for record in tqdm(SeqIO.parse(gb_file, "genbank")):
        seq_num += 1
        SeqIO.write(record, os.path.join(out_dir, "seq_%d.gb" % seq_num), "genbank")
    return seq_num

def gb2fasta(gb_file, out_file)->int:
    """
    convert genbank file to fasta file
    :param gb_file: genbank file
    :param out_file: output fasta file
    :return: sequence count
    """
    seq_num = 0
    with io.open(out_file, "w") as f:
        for record in tqdm(SeqIO.parse(gb_file, "genbank")):
            seq_num += 1
            seq_ana = f">{record.id} {record.description}"
            f.write(formatfasta(seq_ana, str(record.seq)))
    return seq_num

def splitgb_CDS(gb_file, out_file)->int:
    """
    split all CDS to single fasta file from genbank file 
    :param gb_file: genbank file
    :param out_file: output fasta file
    :return: sequence count
    """
    cds_num = 0
    with io.open(out_file, "w") as f:
        for record in tqdm(SeqIO.parse(gb_file, "genbank")):
            for feature in record.features:
                if feature.type == "CDS":
                    cds_num += 1
                    id = feature.qualifiers['protein_id'][0] if 'protein_id' in feature.qualifiers else '<unknown protein_id>'
                    seq_ana = f">lcl|{record.id}_cds_{id}_{cds_num}_{feature.qualifiers['product'][0]}"
                    f.write(formatfasta(seq_ana, str(feature.location.extract(record).seq)))
    return cds_num

def splitgb_CDS_translation(gb_file, out_file)->int:
    """
    split all CDS translation to single fasta file from genbank file
    :param gb_file: genbank file
    :param out_file: output fasta file
    :return: sequence count
    """
    cds_num = 0
    with io.open(out_file, "w") as f:
        for record in tqdm(SeqIO.parse(gb_file, "genbank")):
            for feature in record.features:
                if feature.type == "CDS":
                    cds_num += 1
                    id = feature.qualifiers['protein_id'][0] if 'protein_id' in feature.qualifiers else '<unknown protein_id>'
                    seq_ana = f">lcl|{record.id}_cds_{id}_{cds_num}_{feature.qualifiers['product'][0]}"
                    f.write(formatfasta(seq_ana, str(feature.qualifiers['translation'][0])))
    return cds_num

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Split fasta or gb file into multiple fasta files")
    parser.add_argument("-i", "--input", required=True, help="input file")
    parser.add_argument("-o", "--output", required=True, help="output file or directory")
    parser.add_argument("-t", "--type", required=True, help="input file type, fasta or genbank")
    parser.add_argument("-c", "--cds", action="store_true", help="split all CDS to single fasta file from genbank file")
    parser.add_argument("-s", "--single", action="store_true", help="convert genbank file to single fasta file")
    parser.add_argument("-ct", "--cds_translation", action="store_true", help="split all CDS translation to single fasta file from genbank file")
    args = parser.parse_args()

    if args.type == "fasta":
        seq_num = splitfasta(args.input, args.output)
    elif args.type == "genbank":
        if args.cds:
            seq_num = splitgb_CDS(args.input, args.output)
        elif args.cds_translation:
            seq_num = splitgb_CDS_translation(args.input, args.output)
        elif args.single:
            seq_num = gb2fasta(args.input, args.output)
        else:
            seq_num = splitgb(args.input, args.output)
    else:
        raise ValueError("input file type must be fasta or genbank")

    print(f"split {args.type} file {args.input} into {seq_num} files in directory {args.output}")