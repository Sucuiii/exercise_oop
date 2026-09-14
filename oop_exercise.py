import csv
# Task 1: 基礎類別 GenomicFeature
class GenomicFeature:
    def __init__(self, chromosome: str, start: int, end: int, strand: str):
        if not isinstance(chromosome, str):
            raise ValueError("chromosome must be a string")
        if not (isinstance(start, int) and isinstance(end, int)):
            raise ValueError("start and end must be integers")
        if start <= 0 or end <= 0:
            raise ValueError("start and end must be positive integers")
        if start > end:
            raise ValueError(f"start ({start}) must be <= end ({end})")
        if strand not in ("+", "-"):
            raise ValueError("strand must be '+' or '-'")

        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand

    def length(self) -> int:
        return self.end - self.start + 1

    def overlaps(self, other) -> bool:
        if self.chromosome != other.chromosome:
            return False
        return self.start <= other.end and other.start <= self.end

    def describe(self) -> str:
        return f"{type(self).__name__} {self.chromosome}:{self.start}-{self.end}({self.strand})"


# Task 2 & 3: Exon 子類別
class Exon(GenomicFeature):
    def __init__(self, chromosome: str, start: int, end: int, strand: str, exon_number: int):
        super().__init__(chromosome, start, end, strand)
        self.exon_number = int(exon_number)

    def describe(self) -> str:
        return f"Exon {self.chromosome}:{self.start}-{self.end}({self.strand}) exon #{self.exon_number}"


# Task 3: Gene 子類別
class Gene(GenomicFeature):
    def __init__(self, chromosome: str, start: int, end: int, strand: str, name: str):
        super().__init__(chromosome, start, end, strand)
        self.name = name
        self.exons = []

    def add_exon(self, exon: Exon):
        self.exons.append(exon)

    def total_exon_length(self) -> int:
        return sum(exon.length() for exon in self.exons)

    def describe(self) -> str:
        return f"Gene {self.name} {self.chromosome}:{self.start}-{self.end}({self.strand}), {len(self.exons)} exon(s)"


# Task 3: Variant 子類別
class Variant(GenomicFeature):
    def __init__(self, chromosome: str, start: int, end: int, strand: str, ref_allele: str, alt_allele: str):
        super().__init__(chromosome, start, end, strand)
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele

    def variant_type(self) -> str:
        len_ref = len(self.ref_allele)
        len_alt = len(self.alt_allele)

        if len_ref == 1 and len_alt == 1:
            return "SNP"
        elif len_alt > len_ref:
            return "insertion"
        elif len_alt < len_ref:
            return "deletion"
        else:
            return "MNV"

    def describe(self) -> str:
        v_type = self.variant_type()
        return f"Variant {self.chromosome}:{self.start}-{self.end}({self.strand}) {self.ref_allele}>{self.alt_allele} ({v_type})"


# Task 3: 讀取資料與產生報告
def main():
    genes_dict = {}   # 用來依名稱配對 Exon
    genes = []        # 存放所有 Gene 物件
    variants = []     # 存放所有 Variant 物件

    # 1. 讀取 oop_data.tsv
    with open("oop_data.tsv", "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue

            rec_type, chrom, start_str, end_str, strand = row[0], row[1], row[2], row[3], row[4]
            start = int(start_str)
            end = int(end_str)
            field_a = row[5] if len(row) > 5 else ""
            field_b = row[6] if len(row) > 6 else ""

            if rec_type == "gene":
                gene = Gene(chrom, start, end, strand, name=field_a)
                genes.append(gene)
                genes_dict[field_a] = gene
            elif rec_type == "exon":
                exon = Exon(chrom, start, end, strand, exon_number=int(field_b))
                parent_gene_name = field_a
                if parent_gene_name in genes_dict:
                    genes_dict[parent_gene_name].add_exon(exon)
            elif rec_type == "variant":
                variant = Variant(chrom, start, end, strand, ref_allele=field_a, alt_allele=field_b)
                variants.append(variant)

    # 2. 多型性印出報告 (Polymorphic describe)
    print("=== Features Report ===")
    combined_features = genes + variants
    for feature in combined_features:
        if isinstance(feature, Gene):
            print(f"{feature.describe()} | Total exon length: {feature.total_exon_length()} bp")
        else:
            print(feature.describe())

    # 3. 檢查每個 Variant 落在哪個區域（區分：外側 intergenic、exon 內、intron/UTR 內）
    print("\n=== Variant Locations ===")
    for var in variants:
        overlap_found = False
        for gene in genes:
            if var.overlaps(gene):
                overlap_found = True
                # 檢查是否落在該基因的 exon 內
                in_exon = any(var.overlaps(exon) for exon in gene.exons)
                if in_exon:
                    print(f"{var.describe()} -> Inside {gene.name} (inside exon)")
                else:
                    print(f"{var.describe()} -> Inside {gene.name} (intron / UTR)")
        
        if not overlap_found:
            print(f"{var.describe()} -> intergenic")


if __name__ == "__main__":
    main()