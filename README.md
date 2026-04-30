# E. coli Ampicillin Resistance — Lab Journal
**Курс:** Анализ геномных данных  
**Проект:** #2 — Что вызывает устойчивость к антибиотикам?    

---


### Подготовка рабочей директории

```bash
cd /tmp
mkdir -p project1/raw_data
cd project1/raw_data
```

### Скачивание данных

Референсный геном E. coli K-12 MG1655 (NCBI FTP):

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz
```

Риды резистентного штамма (Google Drive через gdown):

```bash
pip install gdown
gdown --folder "https://drive.google.com/drive/folders/1WhO_qGA2RuEQtTeeC8vh51gO_Dwnh9Km" -O .
unzip amp_res_1.fastq.zip
unzip amp_res_2.fastq.zip
rm *.zip *.gz.1
```

Итог:
-rw-rw-r-- amp_res_1.fastq   122M
-rw-rw-r-- amp_res_2.fastq   122M
-rw-rw-r-- GCF_000005845.2_ASM584v2_genomic.fna.gz   1.4M
-rw-rw-r-- GCF_000005845.2_ASM584v2_genomic.gbff.gz  3.3M
-rw-rw-r-- GCF_000005845.2_ASM584v2_genomic.gff.gz   398K
---

### Ручная инспекция данных

```bash
head -20 amp_res_1.fastq
head -20 amp_res_2.fastq
zcat GCF_000005845.2_ASM584v2_genomic.fna.gz | head -20
wc -l amp_res_1.fastq   # 1823504 → 455876 ридов
wc -l amp_res_2.fastq   # 1823504 → 455876 ридов
```

Формат FASTQ корректный: 4 строки на рид (@header, sequence, +, quality).  
Референс в формате FASTA: одна строка заголовка (>), последовательность разбита по 80 символов.  
Кодирование: Sanger / Illumina 1.9 (Phred33).

---

### Контроль качества — FastQC

```bash
conda config --add pkgs_dirs /tmp/conda_pkgs
conda config --add envs_dirs /tmp/conda_envs
mamba install -c bioconda fastqc -y

mkdir /tmp/project1/fastqc_results
fastqc raw_data/amp_res_1.fastq raw_data/amp_res_2.fastq -o fastqc_results/
```

Результаты FastQC (до трима):

| Модуль | amp_res_1 | amp_res_2 |
|--------|-----------|-----------|
| Per base sequence quality | 🔴 | 🔴 |
| Per sequence quality scores | ✅ | ✅ |
| Per base sequence content | 🟡 | 🟡 |
| Per sequence GC content | 🟡 | 🟡 |

Падение качества на позициях 80–101, нестабильный нуклеотидный состав на первых ~15 позициях — типичные артефакты Illumina.

---

### Фильтрация ридов — Trimmomatic

```bash
mamba install -c bioconda trimmomatic -y

trimmomatic PE \
  raw_data/amp_res_1.fastq \
  raw_data/amp_res_2.fastq \
  raw_data/amp_res_1P.fastq \
  raw_data/amp_res_1U.fastq \
  raw_data/amp_res_2P.fastq \
  raw_data/amp_res_2U.fastq \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20 -phred33
```

Результат:
Input Read Pairs: 455876
Both Surviving: 446259 (97.89%)
Forward Only: 9216 (2.02%)
Reverse Only: 273 (0.06%)
Dropped: 128 (0.03%)
Проверка через wc -l:
```bash
wc -l raw_data/amp_res_1P.fastq  # 1785036 → 446259 ридов ✅
wc -l raw_data/amp_res_2P.fastq  # 1785036 → 446259 ридов ✅
```

FastQC после трима — флаг per base quality сменился с 🔴 на ✅ для обоих файлов.

---

### Выравнивание на референс — BWA-MEM

```bash
mamba install -c bioconda bwa -y

# Распаковка и индексирование референса
gunzip -k raw_data/GCF_000005845.2_ASM584v2_genomic.fna.gz
bwa index raw_data/GCF_000005845.2_ASM584v2_genomic.fna
# Real time: 3.010 sec

# Выравнивание
cd /tmp/project1
bwa mem raw_data/GCF_000005845.2_ASM584v2_genomic.fna \
  raw_data/amp_res_1P.fastq \
  raw_data/amp_res_2P.fastq > alignment.sam
# Real time: 56.947 sec
```

---

### Обработка BAM — SAMtools

```bash
mamba install -c bioconda samtools -y

samtools view -S -b alignment.sam > alignment.bam
samtools flagstat alignment.bam
```

Результат flagstat:
892518 + 0 paired in sequencing
891391 + 0 primary mapped (99.87%)
888554 + 0 properly paired (99.56%)
979 + 0 singletons (0.11%)

```bash
samtools sort alignment.bam -o alignment_sorted.bam
samtools index alignment_sorted.bam
rm alignment.sam  # удалён для экономии места
```

---

### Вариантный коллинг — VarScan2

```bash
mamba install -c bioconda varscan -y

samtools mpileup -f raw_data/GCF_000005845.2_ASM584v2_genomic.fna \
  alignment_sorted.bam > my.mpileup
# 1 samples in 1 input files

varscan mpileup2snp my.mpileup \
  --min-var-freq 0.5 --variants --output-vcf 1 > VarScan_results.vcf
```

Результат:
4641343 bases in pileup file
9 variant positions (6 SNP, 3 indel)
6 variant positions reported (6 SNP, 0 indel)

---

### Аннотация вариантов — SnpEff

```bash
mamba install -c bioconda snpeff -y

echo "k12.genome : ecoli_K12" > snpEff.config
mkdir -p data/k12
gunzip -k raw_data/GCF_000005845.2_ASM584v2_genomic.gbff.gz
cp raw_data/GCF_000005845.2_ASM584v2_genomic.gbff data/k12/genes.gbk

snpEff build -genbank -v k12 \
  -dataDir /tmp/project1/data \
  -config /tmp/project1/snpEff.config

snpEff ann k12 VarScan_results.vcf \
  -dataDir /tmp/project1/data \
  -config /tmp/project1/snpEff.config > VarScan_results_annotated.vcf
```

### Итоговые мутации

| Позиция | Ген | Тип | Замена а/к |
|---------|-----|-----|------------|
| 93 043 | *ftsI* | Миссенс | Ala544Gly |
| 482 698 | *acrB* | Миссенс | Gln569Leu |
| 852 762 | *rybA* | Интрагенная (нкРНК) | — |
| 1 905 761 | *mntP* | Миссенс | Gly25Asp |
| 3 535 147 | *envZ* | Миссенс | Val241Gly |
| 4 390 754 | *rsgA* | Синонимичная | Ala252Ala |
