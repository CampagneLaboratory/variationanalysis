#This readme demonstrates the steps used to train and test a variation model

#Step 1: create a working directory and download the jars and data.

mkdir ~/variations
cd ~/variations

wget http://gobyweb.apps.campagnelab.org/data/DLSV/model-training-1.0.2-snapshot-bin.jar
wget http://gobyweb.apps.campagnelab.org/data/DLSV/goby.jar

wget http://gobyweb.apps.campagnelab.org/data/DLSV/JGSUBWW-ZXLHITM-pickrell-NA19239_argonne-all-files.zip
wget http://gobyweb.apps.campagnelab.org/data/DLSV/VEYNFJY-HRXABYP-pickrell-NA19239_yale-all-files.zip
unzip JGSUBWW-ZXLHITM-pickrell-NA19239_argonne-all-files.zip -d ./somatic
unzip VEYNFJY-HRXABYP-pickrell-NA19239_yale-all-files.zip -d ./germline

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
java -jar goby.jar -m build-sequence-cache human_g1k_v37.fasta.gz

#Step 2: generate an .sbi file with goby
# the .sbi file combines samples into a serialized version of each position,
# which will be used to generate training/testing examples
java -jar goby.jar -m discover-sequence-variants germline/VEYNFJY-HRXABYP-pickrell-NA19239_yale.header somatic/JGSUBWW-ZXLHITM-pickrell-NA19239_argonne.header --format SEQUENCE_BASE_INFORMATION -o ./fullset --genome human_g1k_v37


#Step 2:
#Plant mutations into the set
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.intermediaries.Mutator2 fullset.sbi mutset.sbi

#split the jar file into test, training, and validation here
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.intermediaries.SplitFile -i mutset.sbi -f 0.8 -f 0.1 -f 0.1 -o "set_" -s train -s val -s test

#Step 3: train model with training set
# to save time, we will only traing for 10 epochs (10 full iterations over the whole training set)
# but in practice, you might leave out max-epochs and allow early stopping
# to interrupt training when validation set performance starts detiorating.
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.learning.TrainSomaticModel -t set_train.sbi -v set_val.sbi --max-epochs 10


#Step 4:
# Let's test our model on the test set.
# The AUC performance on the test set will be printed to standard out
# The output of this test is stored in tests/[timestamp]/latest-test.tsv .
# It can be used in conjuction with our MetaR software to generate ROC and Reliability curves.
mkdir tests
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.learning.PredictMutations -i set_test.sbi --long-report -m ./models/[timestamp (eg 1473798870063)]


#Step 5:
# In practice, you will want to use the model on data with no planted mutations,
# that way you can find actual variants. You can use Goby to do this.
# We'll download two pickrell 2010 datasets and look for variants
wget http://dl.dropbox.com/u/357497/RRHCQKJ-discover-sequence-variants-demo-files.zip
unzip RRHCQKJ-discover-sequence-variants-demo-files -d ./practice
#and create a coviarates file to relate two samples to one another
echo -e "sample-id\tpatient-id\tgender\ttype\tkind-of-sample\ttissue\tparents\nPVIZVB-pickrellNA18486_argonne\tP1\tMale\tPatient\tGermline\tBlood\tN/A\nPJCBGUJ-pickrellNA18486_yale\tP1\tMaletatient\tSomatic\tBlood\tN/A" > ./practice/covariates.txt
#Now produce a tsv file with all positions obtaining a variant score (note: not variant probability) for each position
#(TODO: goby3 does not seem to accept custom model path?)
java -jar goby.jar -m discover-sequence-variants ./practice/ZPVIZVB-pickrellNA18486_argonne.header ./practice/PJCBGUJ-pickrellNA18486_yale.header --format SOMATIC_VARIATIONS -o ./practice/NA18486_variants.vcf --genome human_g1k_v37 --covariates ./practice/covariates.txt -X SomaticVariationOutputFormat:model-path="./models/[timestamp]




