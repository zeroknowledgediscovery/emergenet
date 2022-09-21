
for i in `ls *dot`
do

#cat $i | grep -v Prob | sed 's/^Err[^,]*,/",/g' > proc"$i"
cat $i | sed 's/^Err[^,]*,/",/g' > proc"$i"

j=proc"$i"

sed -i 's/Prob:/frac:/g' $j

sed -i 's/\(label=\("\)\)\([0-9]\)/\1index:\n\3/g' $j

dot -Tpdf proc"$i" -o ${j/dot/pdf}
done

