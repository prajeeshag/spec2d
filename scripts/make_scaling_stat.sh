

days=20

grep "Atmos " AMFI_*/stdout_AMFI_* | sed 's/\// /g' > Atmos_time
grep "Atmos_init1step" AMFI_*/stdout_AMFI_* | sed 's/\// /g' > Atmos_time1

npes=($(cat Atmos_time | awk -F "_" '{print $4}'))

xpe=($(cat Atmos_time | sed 's/_/ /g' | awk '{print $5}' | awk -F "x" '{print $1}' ))
ype=($(cat Atmos_time | sed 's/_/ /g' | awk '{print $5}' | awk -F "x" '{print $2}' ))

time1=($(cat Atmos_time | awk '{print $3}'))
time2=($(cat Atmos_time1 | awk '{print $3}'))

echo npes, xpe, ype, time1, time2, time, nyrs_per_day > scaling_stat.csv
for i in "${!npes[@]}"; do
	time=$(echo "${time1[$i]}-${time2[$i]}" | bc)
	nyrs=$(echo "(24.0*3600.0)/(365.0*${time}/${days})" | bc)
	echo ${npes[$i]}, ${xpe[$i]}, ${ype[$i]}, ${time1[$i]}, ${time2[$i]}, $time, $nyrs >> scaling_stat.csv
done
