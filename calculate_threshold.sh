file=$1
java -Xmx5G -jar Aracne.jar -e $file -o . --tfs TFs_list.txt --pvalue 1E-8 --seed 1 --calculateThreshold

