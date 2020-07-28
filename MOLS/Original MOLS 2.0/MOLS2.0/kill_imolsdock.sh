#shell script to kill imolsdock job
for x in `cat imolsdock.pid`
	do
	kill -15 $x
	done
