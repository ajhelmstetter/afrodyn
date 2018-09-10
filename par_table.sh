
#!/bin/bash

find . -name genes_with_paralog_warnings.txt >> par_list.txt

rm para_table.txt
touch para_table.txt

while read i
do
cat $i >> para_table.txt
echo $i >> para_table.txt

echo -e "\n" >> para_table.txt
done < par_list.txt

sed -i '/^$/d' para_table.txt
