# split
A tool to split gb or fasta files downloaded in bulk from NCBI.

## 🎯 Requirement
- Python>=3.10
- Biopython
## 🌟 Usage
You can get detailed help with `-h`.

Here are some examples.

**FASTA to FASTAs**
```
python split.py -i seq.fasta -o result -t fasta
```

**GB to GBS**
```
python split.py -i seq.gb -o result -t genbank
```

**GB to FASTAs**
```
python split.py -i seq.gb -o result.fasta -t genbank -s
python split.py -i result.fasta -o result -t fasta
```

**Split CDS from GB**
```
python split.py -i seq.gb -out result -t genbank -c
```