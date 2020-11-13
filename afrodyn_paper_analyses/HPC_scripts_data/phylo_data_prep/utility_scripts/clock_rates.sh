#!/bin/bash

#cp ebsp_anni.xml test.xml

i=0

while read line
do
	echo "$line"
    if [ "$i" -eq "0" ]
    then
        sed -i -e "s|aligned\" estimate=\"false\" name=\"clock.rate\">1.0<|aligned\" estimate=\"false\" name=\"clock.rate\">"$line"<|g" $1
        i=$((i+1))
        echo "$i"
    else
    	sed -i -e "s|aligned"$i"\" estimate=\"false\" name=\"clock.rate\">1.0<|aligned"$i"\" estimate=\"false\" name=\"clock.rate\">"$line"<|g" $1
    	i=$((i+1))
    	echo "aligned$i"
    fi
done < $2
