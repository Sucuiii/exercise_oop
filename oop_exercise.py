class GenomicFeature:
    def __init__(self, chromosome: str, start: int, end: int, strand: str):
        # 1. 驗證 chromosome 必須是字串
        if not isinstance(chromosome, str):
            raise ValueError("chromosome must be a string")

        # 2. 驗證 start 和 end 必須是正整數（1-based），且 start <= end
        if not (isinstance(start, int) and isinstance(end, int)):
            raise ValueError("start and end must be integers")
        if start <= 0 or end <= 0:
            raise ValueError("start and end must be positive integers")
        if start > end:
            raise ValueError(f"start ({start}) must be <= end ({end})")

        # 3. 驗證 strand 必須是 '+' 或 '-'
        if strand not in ("+", "-"):
            raise ValueError("strand must be '+' or '-'")

        # 4. 將屬性存入物件
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand

    def length(self) -> int:
        # 計算特徵涵蓋的 base 數（end - start + 1）
        return self.end - self.start + 1

    def overlaps(self, other) -> bool:
        # 同一條染色體才可能重疊
        if self.chromosome != other.chromosome:
            return False
        # 閉區間重疊判定：A.start <= B.end 且 B.start <= A.end
        return self.start <= other.end and other.start <= self.end

    def describe(self) -> str:
        # type(self).__name__ 會動態抓取目前的類別名稱（GenomicFeature）
        # 後面 Task 2 繼承給 Exon 時就會自動變成 Exon
        cls_name = type(self).__name__
        return f"{cls_name} {self.chromosome}:{self.start}-{self.end}({self.strand})"


if __name__ == "__main__":
    # 這是 README 和投影片給你的測試範例
    a = GenomicFeature("chr1", 1000, 5000, "+")
    b = GenomicFeature("chr1", 4800, 6000, "+")
    c = GenomicFeature("chr2", 1000, 5000, "+")

    print(a.describe())     # 印出: GenomicFeature chr1:1000-5000(+)
    print(a.length())       # 印出: 4001
    print(a.overlaps(b))    # 印出: True
    print(a.overlaps(c))    # 印出: False

    # 測試異常輸入是否會正確拋出 ValueError
    try:
        GenomicFeature("chr1", 5000, 1000, "+")
    except ValueError as e:
        print(f"Correctly caught error: {e}")