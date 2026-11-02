#!/usr/bin/env python3
"""
make_test_data.py -- synthetic end-to-end test dataset (no internet, no real data).

Creates in tests/data/:
  taxdump/        mini NCBI-style taxonomy (fake taxids >= 900000, real ones for human/bacteria)
  refs.fasta + refs.taxid.map      custom ITS "database" (6 plant species + decoys)
  core_nt_mock.fasta + .taxid.map  stand-in for core_nt (refs + human, E. coli, a fungus)
  human_ref.fasta                  tiny stand-in "human genome"
  raw/<sample>_L001_R{1,2}_001.fastq.gz   simulated Illumina 2x150 pairs

ITSx needs real conserved boundaries, so the ITS constructs use the consensus of
ITSx's own boundary HMMs (hmmemit) around random ITS1/ITS2.

Samples and the expected outcome:
  T_A  Alpha one 60% + Beta one 40% (+5% human, +3% E. coli reads)
       -> genus Alpha (Alpha two in the DB differs by 1 SNP per spacer: a TIE, so genus, by design),
          species Beta one (Beta two is 5% different: resolved), human row, E. coli via core_nt control
  T_D  Alpha two only (99.6% identical to Alpha one in the DB)
       -> TIE -> genus Alpha (not a false species)
  T_E  Alpha four = NOT in the database, ~6% diverged from Alpha one/two
       -> identity cap -> genus Alpha (not forced onto its nearest relative)
  T_B  tiny sample (assembly fails): mostly human + a few Alpha one ITS reads
       -> read-level path: human row + genus/species from reads
  T_C  assembles but contains NO ITS (fungal-like random sequence)
       -> fallback contigs, no custom-DB hit, core_nt control names it a fungus
"""
import gzip, os, random, shutil, subprocess, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(ROOT, "data")
random.seed(7)
B = "ACGT"


def rnd(n):
    return "".join(random.choice(B) for _ in range(n))


def mutate(s, n_sub=0, frac=0.0):
    s = list(s)
    k = n_sub if n_sub else int(len(s) * frac)
    for i in random.sample(range(len(s)), k):
        s[i] = random.choice([b for b in B if b != s[i]])
    return "".join(s)


def rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def hmm_consensus():
    hmm = None
    for p in os.environ["PATH"].split(":"):
        cand = os.path.join(p, "ITSx_db", "HMMs", "T.hmm")
        if os.path.exists(cand):
            hmm = cand
            break
    if hmm is None:
        sys.exit("ITSx_db/HMMs/T.hmm not found next to ITSx in PATH")
    out = subprocess.run(["hmmemit", "-c", hmm], capture_output=True, text=True, check=True).stdout
    cons, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].split("-consensus")[0]
            cons[name] = ""
        else:
            cons[name] += line.strip()
    pick = lambda frag: next(v for k, v in cons.items() if frag in k).upper().replace("U", "T")
    return (pick("1_SSU_Tracheophyta_1_SSU_end_one_long"),
            pick("2_5.8_Tracheophyta_58S_start_onea_long"),
            pick("3_End_Tracheophyta_58S_end_one_long"),
            pick("4_LSU_Tracheophyta_4_LSU_start_one_long"))


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for h, s in recs:
            fh.write(f">{h}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


def fq_gz(path, reads):
    with gzip.open(path, "wt") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def simulate_pairs(template, n_pairs, tag, rl=150, err=0.003):
    r1, r2 = [], []
    L = len(template)
    for i in range(n_pairs):
        ins = max(int(random.gauss(320, 90)), rl + 10)
        ins = min(ins, L)
        start = random.randint(0, L - ins)
        frag = template[start:start + ins]
        if random.random() < 0.5:
            frag = rc(frag)
        a, b = frag[:rl], rc(frag)[:rl]
        a, b = mutate(a, frac=err) if err else a, mutate(b, frac=err) if err else b
        r1.append((f"{tag}.{i} 1:N:0:1", a))
        r2.append((f"{tag}.{i} 2:N:0:1", b))
    return r1, r2


def main():
    shutil.rmtree(D, ignore_errors=True)
    os.makedirs(f"{D}/taxdump"); os.makedirs(f"{D}/raw"); os.makedirs(f"{D}/human")
    ssu_end, s58_start, s58_end, lsu_start = hmm_consensus()
    mid58 = rnd(70)

    def its(its1, its2):
        return its1 + s58_start + mid58 + s58_end + its2

    # --- ITS regions of the "database" species ------------------------------------
    A_its1, A_its2 = rnd(230), rnd(220)
    alpha1 = its(A_its1, A_its2)
    alpha2 = its(mutate(A_its1, n_sub=1), mutate(A_its2, n_sub=1))           # ~99.6% id to alpha1
    alpha3 = its(mutate(A_its1, frac=0.04), mutate(A_its2, frac=0.04))
    alpha4 = its(mutate(A_its1, frac=0.07), mutate(A_its2, frac=0.07))       # NOT in the db (novel)
    B_its1, B_its2 = rnd(240), rnd(210)
    beta1 = its(B_its1, B_its2)
    beta2 = its(mutate(B_its1, frac=0.05), mutate(B_its2, frac=0.05))
    gamma1 = its(mutate(B_its1, frac=0.14), mutate(B_its2, frac=0.14))       # other genus, same family

    refs = [  # id, taxid, seq
        ("REF001.1_full", 900003, alpha1), ("REF002.1_full", 900004, alpha2),
        ("REF003.1_full", 900005, alpha3), ("REF004.1_full", 900031, alpha1),   # subspecies of alpha one
        ("REF005.1_full", 900013, beta1), ("REF006.1_full", 900014, beta2),
        ("REF007.1_full", 900021, gamma1),
    ]
    human_its = rnd(600)
    ecoli_seq = rnd(1500)
    fungus_seq = rnd(2000)                                     # no ITS HMM signal
    write_fasta(f"{D}/refs.fasta", [(i, s) for i, _, s in refs] + [("DECOY1.1_full", human_its)])
    with open(f"{D}/refs.taxid.map", "w") as fh:
        for i, t, _ in refs:
            fh.write(f"{i}\t{t}\n")
        fh.write("DECOY1.1_full\t9606\n")
    write_fasta(f"{D}/core_nt_mock.fasta",
                [(i, s) for i, _, s in refs] + [("DECOY1.1_full", human_its), ("ECOLI1.1", ecoli_seq),
                                                 ("FUNG1.1", fungus_seq)])
    with open(f"{D}/core_nt_mock.taxid.map", "w") as fh:
        for i, t, _ in refs:
            fh.write(f"{i}\t{t}\n")
        fh.write(f"DECOY1.1_full\t9606\nECOLI1.1\t562\nFUNG1.1\t900050\n")

    # --- mini taxonomy ---------------------------------------------------------------------
    nodes = [  # taxid, parent, rank, name
        (1, 1, "no rank", "root"), (131567, 1, "no rank", "cellular organisms"),
        (2, 131567, "superkingdom", "Bacteria"), (2759, 131567, "superkingdom", "Eukaryota"),
        (33090, 2759, "kingdom", "Viridiplantae"), (4751, 2759, "kingdom", "Fungi"),
        (33208, 2759, "kingdom", "Metazoa"),
        (9604, 33208, "family", "Hominidae"), (9605, 9604, "genus", "Homo"), (9606, 9605, "species", "Homo sapiens"),
        (543, 2, "family", "Enterobacteriaceae"), (561, 543, "genus", "Escherichia"), (562, 561, "species", "Escherichia coli"),
        (900001, 33090, "family", "Testaceae"), (900002, 900001, "genus", "Alpha"),
        (900003, 900002, "species", "Alpha one"), (900004, 900002, "species", "Alpha two"),
        (900005, 900002, "species", "Alpha three"), (900031, 900003, "subspecies", "Alpha one subsp. x"),
        (900011, 33090, "family", "Betaceae"), (900012, 900011, "genus", "Beta"),
        (900013, 900012, "species", "Beta one"), (900014, 900012, "species", "Beta two"),
        (900020, 900011, "genus", "Gamma"), (900021, 900020, "species", "Gamma one"),
        (900049, 4751, "genus", "Fungusia"), (900050, 900049, "species", "Fungusia sp."),
    ]
    with open(f"{D}/taxdump/nodes.dmp", "w") as fn, open(f"{D}/taxdump/names.dmp", "w") as fm:
        for t, p, r, n in nodes:
            fn.write(f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\n")
            fm.write(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n")

    # --- human stand-in reference ------------------------------------------------------------
    human_ref = rnd(40000)
    write_fasta(f"{D}/human/human_ref.fasta", [("chrTEST", human_ref)])

    # --- samples ------------------------------------------------------------------------------
    def cistron(its_seq):
        return rnd(700) + ssu_end + its_seq + lsu_start + rnd(700)

    def human_pairs(n, tag):
        return simulate_pairs(human_ref, n, tag)

    def build(sample, parts):
        r1, r2 = [], []
        for tag, tmpl, n in parts:
            a, b = simulate_pairs(tmpl, n, f"{sample}.{tag}")
            r1 += a; r2 += b
        idx = list(range(len(r1)))
        random.shuffle(idx)
        r1 = [r1[i] for i in idx]; r2 = [r2[i] for i in idx]
        # unique names per sample
        fq_gz(f"{D}/raw/{sample}_L001_R1_001.fastq.gz", r1)
        fq_gz(f"{D}/raw/{sample}_L001_R2_001.fastq.gz", r2)

    ssu_a1 = cistron(alpha1)
    ssu_b1 = cistron(beta1)
    build("T_A", [("alpha1", ssu_a1, 3600), ("beta1", ssu_b1, 2400), ("human", human_ref, 300), ("ecoli", ecoli_seq, 180)])
    build("T_D", [("alpha2", cistron(alpha2), 5000)])
    build("T_E", [("alpha4", cistron(alpha4), 5000)])
    build("T_B", [("human", human_ref, 130), ("alpha1", ssu_a1, 20)])
    build("T_C", [("fungus", fungus_seq, 4000)])
    print("synthetic dataset written to", D)


if __name__ == "__main__":
    main()
