#scripts QIIME2. 2022.2
#Dr. Juan Alfredo Hernández García
#Instituto Politécnico Nacional
#Escuela Nacional de Ciencias Biológicas
#Dpartamento de Microbiología

`conda activate qiime2-2022.2.

echo 'alias qiime2="source activate qiime2-2022.2"' >> $HOME/.bashrc
echo 'alias picrust2="source activate picrust2"' >> $HOME/.bashrc
source $HOME/.bashrc`

#actualizamos conda
`conda update conda`

#enlistamos los comandos para conocer los ambientes de conda
`conda env list`

#Actualizar todos los paquetes de conda (base)
`conda update conda --all`

#Eliminar temporales y actualizar todo
`conda clean --al
conda update --all`
_______________________________________________________________________________
# **Preprocesamiento de datos metabarcoding 16S rRNA**
#1.0 Obtención de datos
## Tres posibles escenarios:
MiSeq, un archivo por muestra
`sampleA_S1_R1_001.fastq.gz
sampleA_S1_R2_001.fastq.gz`

find -maxdepth 1 -name "*fastq.gz" -type f | rename 's/_001//'

##NextSeq, cuatro archivos por muestra, L001-L004
```
sampleA_S1_L001_R1_001.fastq.gz
sampleA_S1_L001_R2_001.fastq.gz
sampleA_S1_L002_R1_001.fastq.gz
sampleA_S1_L002_R2_001.fastq.gz
sampleA_S1_L003_R1_001.fastq.gz
sampleA_S1_L003_R2_001.fastq.gz
sampleA_S1_L004_R1_001.fastq.gz
sampleA_S1_L004_R2_001.fastq.gz
```

##Podemos concatenar los archivos unos por uno :(
`cat sampleA_S1_L00?_R1* > sampleA_S1_L001_R1_001.fastq.gz
cat sampleA_S1_L00?_R2* > sampleA_S1_L001_R2_001.fastq.gz`
## O podemos usar iteraciones para concatenar los archivos de las 4 líneas, L001-L004, de Illumina al formato `sampleA_S1_R1.fastq.gz` con el siguiente LOOP

```
while read i
do
cat ${i}*R1*.fastq.gz > ${i}_R1.fastq.gz
cat ${i}*R2*.fastq.gz > ${i}_R2.fastq.gz
done < <(ls *R1*gz | cut -d\_ -f1,2 | sort | uniq )
```

#C) NextSeq, un archivo por muestra, L001
```
sampleA_S1_L001_R1_001.fastq.gz
sampleA_S1_L001_R2_001.fastq.gz
```
______________________________________________________________________________
# Para los pasos siguientes pasos utilizaremos el siguiente set de prueba con lecturas pair-end de Illumina.
``wget https://datos16s.s3.us-east-2.amazonaws.com/atg16s.tar.gz
tar -zxvf atg16s.tar.gz && rm atg16s.tar.gz``

#**Control de calidad**
Cada archivo fastq tiene cuatro lineas:
1. Nombre de la secuencia (header - id del secuenciador, coordenadas del spot, flow-cell,
adaptador, etc.)
2. Secuencia
3. Espaciador (+)
4. Valores de calidad: Q Score - alfanumérico.
Para evaluar la calidad se podría evaluar los archivos de manera individual pero eso consumiría
"mucho tiempo", por lo que es mejor tratarlo como un sólo archivo, ¿cómo lo haríamos?...

Concatenamos todos los archivos .fastq en un solo y ejecutamos el comando `fastqc`.
Creamos un folder para realizar el análisis de QIIME
```
mkdir 01_qc
source activate qc
cat *.fastq.gz > all.fq.gz
fastqc --nogroup -f fastq all.fq.gz
conda deactivate
# NOTA: Borrar el archivo concatenado anteriormente: :
rm all_fastqc.zip all.fq.gz
```

# **IMPORTAR DATOS*
##Es necesario crear un archivo manifest donde se especifique sample-id, absolute-filepath y direction por cada uno de las muestras.
# Colocamos todas nuestras secuencias en el directorio de dataset y ejecutamos el siguiente LOOP

```
find dataset/ -name "*fastq.gz" -type f | rename 's/_L001//; s/_001//'
dir="dataset"
ids=$(ls ${dir}/*gz | cut -d\_ -f1,2 | sort | uniq)
echo "sample-id,absolute-filepath,direction" > manifest.csv
for i in ${ids}
do
name=${i#*/}; name=${name%_*}
echo "${name},\$PWD/${i}_1.fastq.gz,forward" >> manifest.csv
echo "${name},\$PWD/${i}_2.fastq.gz,reverse" >> manifest.csv
done
```

#Se puede utilizar un plugin de Google Docs Sheets para verificar la integridad de nuestro archivo de metadatos: http://keemei.qiime.org/
#QIIME2 utliza artefactos: * Artifacto = fastq + manifest * Tienen la extensión qza, qiime zip artifact
#Si obtenemos algún error de conda en el entorno qiime2 en AWS.
```
mkdir -p 01_qc
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.csv --input-format PairedEndFastqManifestPhred33 \
--output-path 01_qc/01_demux.qza
```
## Visualización
```
qiime demux summarize \
--i-data 01_qc/01_demux.qza \
--o-visualization 01_qc/01_demux.qzv
```
#Los archivos qzv son artefactos de visualización, se puede ver su contenido de manera local con:
`qiime tools view 01_qc/01_demux.qzv #te redirecciona a una página web`

#**Quitar adaptadores**
# 338F ACTCCTACGGGAGGCAGCA \
# 806R GGACTACHVGGGTWTCTAAT
```
qiime cutadapt trim-paired \
--i-demultiplexed-sequences 01_qc/01_demux.qza \
--p-cores "$(nproc)" \
--p-front-f ACTCCTACGGGAGGCAGCA \
--p-front-r GGACTACHVGGGTWTCTAAT \
--o-trimmed-sequences 01_qc/01_demux-trim.qza \
--verbose
```
qiime demux summarize \
--i-data 01_qc/01_demux-trim.qza \
--o-visualization 01_qc/01_demux-trim.qzv

`qiime tools view 01_qc/01_demux-trim.qzv #te redirecciona a una página web`


#**Denoising - Eliminación de ruido**
Los pasos que realiza DADA2 son:
* Filtrar y recortar
* Dereplicar
* Aprendizaje de tazas de error
* Inferencia de la composición de la muestra
* Unión de lecturas paired-end
* Creación de tabla de secuencias
* Eliminación de quimeras

```
qiime dada2 denoise-paired \
--i-demultiplexed-seqs 01_qc/01_demux.qza \
--output-dir 02_dada2.1 \
--p-trim-left-f 9 \
--p-trim-left-r 9 \
--p-trunc-len-f 200 \
--p-trunc-len-r 200 \
--p-n-threads "$(nproc)"
```
#a) Denoising stats
`qiime metadata tabulate \
--m-input-file 02_dada2.1/denoising_stats.qza \
--o-visualization 02_dada2/stats.qzv`
`qiime tools view 02_dada2/stats.qzv #te redirecciona a una página web`
### * el archivo _metadata.tsv_ se descarga de la página web

#b) Table
`qiime feature-table summarize \
--i-table 02_dada2.1/table.qza \
--o-visualization 02_dada2.1/table.qzv \
--m-sample-metadata-file 02_dada2.1/metadata.tsv`

`qiime tools view 02_dada2.1/table.qzv #te redirecciona a una página web`

#c) Representative sequences
`qiime feature-table tabulate-seqs \
--i-data 02_dada2.1/representative_sequences.qza \
--o-visualization 02_dada2.1/rep-seqs.qzv`

`qiime tools view 02_dada2/rep-seqs.qzv #te redirecciona a una página web`

# **Clasificación taxonómica**
#1.- Usaremos greengenes data base
#Pre-fitted sklearn-based taxonomy classifier

`wget https://data.qiime2.org/2020.8/common/gg-13-8-99-515-806-nb-classifier.qza
mkdir 03_taxonomy

qiime feature-classifier classify-sklearn \
--i-reads 02_dada2/representative_sequences.qza \
--i-classifier training_feature_classifiers/classifier.qza \
--o-classification 03_taxonomy/03_dada2-taxonomy-sklearn.qza \
--verbose

qiime metadata tabulate \
--m-input-file 03_taxonomy/03_dada2-taxonomy-sklearn.qza \
--o-visualization 03_taxonomy/03_dada2-taxonomy-sklearn.qzv

qiime taxa barplot \
--i-table 02_dada2/table.qza \
--i-taxonomy 03_taxonomy/03_dada2-taxonomy-sklearn.qza \
--m-metadata-file metadata.tsv \
--o-visualization 03_taxonomy/03_dada2-taxonomy-sklearn-barplots.qzv

# Agrupar tabla por categoría
qiime feature-table group \
--i-table 02_dada2/table.qza \
--p-axis sample \
--m-metadata-file metadata.tsv \
--m-metadata-column Treatment \
--p-mode sum \
--o-grouped-table 03_taxonomy/03_dada2-tax-sk_Treatment.qza

qiime metadata tabulate \
--m-input-file 03_taxonomy/03_dada2-tax-sk_Treatment.qza \
--o-visualization 03_taxonomy/03_dada2-tax-sk_Treatment.qzv`
`qiime tools view 03_taxonomy/03_dada2-tax-sk_Treatment.qzv #te redirecciona a una página web`


`qiime feature-table group \
--i-table 02_dada2/table.qza \
--p-axis sample \
--m-metadata-file metadata.tsv \
--m-metadata-column Gut_fraction \
--p-mode sum \
--o-grouped-table 03_taxonomy/03_dada2-tax-sk_GUT.qza`

`qiime metadata tabulate \
--m-input-file 03_taxonomy/03_dada2-tax-sk_GUT.qza \
--o-visualization 03_taxonomy/03_dada2-tax-sk_GUT.qzv`
`qiime tools view 03_taxonomy/03_dada2-tax-sk_GUT.qzv #te redirecciona a una página web`


#Extraemos los datos de la columna 2 y 3, específicamente del tratamiento
``echo 'id' > treatment.txt; awk '{print $2}' metadata.tsv | \
sort | uniq >> treatment.txt`

`qiime taxa barplot \
--i-table 03_taxonomy/03_dada2-tax-sk_Treatment.qza \
--i-taxonomy 03_taxonomy/03_dada2-taxonomy-sklearn.qza \
--m-metadata-file treatment.txt \
--o-visualization 03_taxonomy/03_dada2-tax-sk-barplots_Treatment.qzv
`qiime tools view 03_taxonomy/03_dada2-tax-sk-barplots_Treatment.qzv #te redirecciona a una página web`

`echo 'id' > treatment.txt; awk '{print $2}' metadata.tsv | \
sort | uniq >> treatment.txt`

`qiime taxa barplot \
--i-table 03_taxonomy/03_dada2-tax-sk_Treatment.qza \
--i-taxonomy 03_taxonomy/03_dada2-taxonomy-sklearn.qza \
--m-metadata-file treatment.txt \
--o-visualization 03_taxonomy/03_dada2-tax-sk-barplots_Treatment.qzv
`qiime tools view 03_taxonomy/03_dada2-tax-sk-barplots_Treatment.qzv #te redirecciona a una página web


#2.- Usaremos **SILVA** data base

#Pre-fitted sklearn-based taxonomy classifier

mkdir 03_taxonomy_SILVA

qiime feature-classifier classify-sklearn \
--i-reads 02_dada2/representative_sequences.qza \
--i-classifier silva-138-99-nb-classifier.qza \
--o-classification 03_taxonomy_SILVA/03_dada2-taxonomy-sklearn.qza \
--verbose

qiime metadata tabulate \
--m-input-file 03_taxonomy_SILVA/03_dada2-taxonomy-sklearn.qza \
--o-visualization 03_taxonomy_SILVA/03_dada2-taxonomy-sklearn.qzv

qiime taxa barplot \
--i-table 02_dada2/table.qza \
--i-taxonomy 03_taxonomy_SILVA/03_dada2-taxonomy-sklearn.qza \
--m-metadata-file metadata.tsv \
--o-visualization 03_taxonomy_SILVA/03_dada2-taxonomy-sklearn-barplots.qzv

# Agrupar tabla por categoría
qiime feature-table group \
--i-table 02_dada2/table.qza \
--p-axis sample \
--m-metadata-file metadata.tsv \
--m-metadata-column Treatment \
--p-mode sum \
--o-grouped-table 03_taxonomy_SILVA/03_dada2-tax-sk_Treatment.qza

qiime metadata tabulate \
--m-input-file 03_taxonomy_SILVA/03_dada2-tax-sk_Treatment.qza \
--o-visualization 03_taxonomy_SILVA/03_dada2-tax-sk_Treatment.qzv
`qiime tools view 03_taxonomy_SILVA/03_dada2-tax-sk_Treatment.qzv #te redirecciona a una página web`

#Extraemos los datos de la columna 2, específicamente del tratamiento
echo 'id' > treatment.txt; awk '{print $2}' metadata.tsv | \
sort | uniq >> treatment.txt

qiime taxa barplot \
--i-table 03_taxonomy/03_dada2-tax-sk_Treatment.qza \
--i-taxonomy 03_taxonomy/03_dada2-taxonomy-sklearn.qza \
--m-metadata-file treatment.txt \
--o-visualization 03_taxonomy_SILVA/03_dada2-tax-sk-barplots_Treatment.qzv
`qiime tools view 03_taxonomy_SILVA/03_dada2-tax-sk-barplots_Treatment.qzv #te redirecciona a una página web`
