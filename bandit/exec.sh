S=0
make && S=1

if [ $S -eq 0 ] ; then
	echo -e '\e[38;5;220m --- make fault --- \e[m'
	exit 0
fi

echo -n -e '\e[38;5;219m こ、これでどうかな ... ?\e[m'
read F

rm -f ^a

if [ -z $F ]; then
	echo -e '\e[38;5;226m --- valgrind --- \e[m'
	rm -f .valgrind.log
	valgrind --leak-check=full --error-exitcode=222 --log-file='.valgrind.log' `grep -v '#' _exec.sh` | tee ^a
	[ `grep "ERROR SUMMARY: 0 errors" .valgrind.log | wc -l` -ne 1 ] \
		&& echo -n -e '\e[38;5;202m --- valgrind error --- \e[m' \
		&& read \
		&& cat .valgrind.log
else
	echo -e '\e[38;5;110m --- raw --- \e[m'
	`grep -v '#' _exec.sh` | tee ^a
fi
