"""Convert pinned GEO supplements to gene-by-sample matrices without web lookups."""

import csv
import gzip
import io
import re
import tarfile
from pathlib import Path

from openpyxl import load_workbook

CACHE = Path(__file__).resolve().parent / "geo_cache"


def write_table(name, header, rows):
    path = CACHE / name
    temporary = path.with_suffix(".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    temporary.replace(path)


def annotation():
    gene_names, tx_to_gene = {}, {}
    with gzip.open(CACHE / "Homo_sapiens.GRCh38.97.gtf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] not in {"gene", "transcript"}:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            gene = attrs["gene_id"]
            if fields[2] == "gene":
                gene_names[gene] = attrs.get("gene_name", "")
            else:
                tx_to_gene[attrs["transcript_id"]] = gene
    write_table("ensembl97_symbols.csv", ["gene_id", "symbol"], sorted(gene_names.items()))
    return gene_names, tx_to_gene


def imatinib():
    path = CACHE / "GSE321679_Normalised_counts_all_samples.xlsx"
    with path.open("rb") as handle:
        book = load_workbook(handle, read_only=True, data_only=True)
        sheet = book["K562-WT- DMSO VS Imatinib"]
        rows = sheet.values
        header = [str(value).strip() for value in next(rows)]
        columns = [i for i, value in enumerate(header) if re.fullmatch(r"WT-(imatinib|DMSO)-[123]", value)]
        if len(columns) != 6:
            raise ValueError(f"Expected six Imatinib sample columns, got {header}")
        write_table("GSE321679_logcpm.csv", ["gene_id", "symbol", *(header[i] for i in columns)],
                    ((row[0], row[1], *(row[i] for i in columns)) for row in rows
                     if row[0] is not None))
        book.close()


def count_matrix(files, name):
    sample_names = list(files)
    first = files[sample_names[0]]
    for sample, records in files.items():
        if records.keys() != first.keys():
            raise ValueError(f"Gene IDs differ across {name} samples: {sample}")
    write_table(name, ["gene_id", "symbol", *sample_names],
                ((gene, first[gene][0], *(data[gene][1] for data in files.values()))
                 for gene in sorted(first)))


def vorinostat():
    files = {}
    with tarfile.open(CACHE / "GSE319257_RAW.tar") as tar:
        for member in tar.getmembers():
            if not re.fullmatch(r"GSM\d+_C\d+_mock_(control|vorinostat)_[1-6]\.xlsx", member.name):
                continue
            with tar.extractfile(member) as handle:
                book = load_workbook(io.BytesIO(handle.read()), read_only=True, data_only=True)
                sheet = book.worksheets[0]
                sheet.reset_dimensions()
                rows = sheet.values
                header = next(rows)
                gene, symbol, counts = (header.index(key) for key in ("Gene ID", "Name", "Unique gene reads"))
                records = {}
                for row in rows:
                    if not row[gene]:
                        continue
                    if row[gene] in records:
                        raise ValueError(f"Duplicate gene ID in {member.name}: {row[gene]}")
                    records[row[gene]] = (row[symbol], row[counts])
                files[member.name.removesuffix(".xlsx")] = records
                book.close()
    count_matrix(files, "GSE319257_counts.csv")


def lps():
    names = {
        "GSM8993640": "donor1_ctrl",
        "GSM8993641": "donor1_treat",
        "GSM8993646": "donor2_ctrl",
        "GSM8993647": "donor2_treat",
    }
    files = {}
    with tarfile.open(CACHE / "GSE297509_RAW.tar") as tar:
        for member in tar.getmembers():
            accession = member.name.split("_")[0]
            if accession not in names:
                continue
            records = {}
            with tar.extractfile(member) as raw, gzip.open(raw, "rt") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    gene = row["Gene ID"]
                    if gene in records:
                        raise ValueError(f"Duplicate gene ID in {member.name}: {gene}")
                    records[gene] = (row["Name"], int(row["unique reads"]))
            files[f"{accession}_{names[accession]}"] = records
    count_matrix(files, "GSE297509_counts.csv")


def tnf(gene_names, tx_to_gene):
    names = {
        "GSM9185037_C156", "GSM9185038_T156", "GSM9185044_C67",
        "GSM9185045_T67", "GSM9185051_C840", "GSM9185052_T840",
    }
    files, qc = {}, []
    with tarfile.open(CACHE / "GSE305796_RAW.tar") as tar:
        for member in tar.getmembers():
            sample = member.name.removesuffix("_quant.sf.gz")
            if sample not in names:
                continue
            gene_tpm = {}
            total_reads = mapped_reads = 0.0
            with tar.extractfile(member) as raw, gzip.open(raw, "rt") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    reads = float(row["NumReads"])
                    total_reads += reads
                    gene = tx_to_gene.get(row["Name"].split(".")[0])
                    if gene:
                        gene_tpm[gene] = gene_tpm.get(gene, 0.0) + float(row["TPM"])
                        mapped_reads += reads
            if mapped_reads / total_reads < 0.95:
                raise ValueError(f"Low Ensembl mapping coverage for {sample}")
            # Scale gene-level TPMs to library size, as in tximport's scaledTPM.
            scale = total_reads / sum(gene_tpm.values())
            files[sample] = {g: (gene_names[g], tpm * scale) for g, tpm in gene_tpm.items()}
            qc.append((sample, total_reads, mapped_reads, mapped_reads / total_reads))
    count_matrix(files, "GSE305796_scaled_tpm_counts.csv")
    write_table("GSE305796_mapping_qc.csv", ["sample", "total_reads", "mapped_reads", "mapped_fraction"], qc)


if __name__ == "__main__":
    symbols, transcripts = annotation()
    imatinib()
    vorinostat()
    lps()
    tnf(symbols, transcripts)
