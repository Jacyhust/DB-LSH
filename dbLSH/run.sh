cmake
make wlsh

:<<EOF
./wlsh audio 1.5 50 5 10 0.1

./wlsh mnist 1.5 50 5 10 0.1

./wlsh cifar 1.5 50 5 10 0.1

./wlsh NUS 1.5 50 5 10 0.1

./wlsh Trevi 1.5 50 5 10 0.1 700

./wlsh deep1m 1.5 50 5 10 0.1 0.2 #optimal

./wlsh gist 1.5 50 5 10 0.1 0.3 #optimal




#./wlsh tiny 1.5 50 5 12 0.15 0.1

#./wlsh tiny 1.5 50 5 12 0.1 0.15

#./wlsh tiny 1.5 50 5 12 0.1 0.15
EOF

./wlsh cifar 1.5 50 5 10 0.1 300
./wlsh cifar 1.5 50 5 10 0.1 250
./wlsh cifar 1.5 50 5 10 0.1 350

./wlsh cifar 1.5 50 5 10 0.15 300
./wlsh cifar 1.5 50 5 10 0.15 250
./wlsh cifar 1.5 50 5 10 0.15 350

./wlsh Trevi 1.5 50 5 10 0.1 700
./wlsh Trevi 1.5 50 5 10 0.1 600
./wlsh Trevi 1.5 50 5 10 0.1 800

./wlsh Trevi 1.5 50 5 10 0.15 700
./wlsh Trevi 1.5 50 5 10 0.15 600 #optimal
./wlsh Trevi 1.5 50 5 10 0.15 800

./wlsh NUS 1.5 50 5 10 0.15 6.5
./wlsh NUS 1.5 50 5 10 0.15 7
./wlsh NUS 1.5 50 5 10 0.15 7.5

./wlsh NUS 1.5 50 5 10 0.1 6.5
./wlsh NUS 1.5 50 5 10 0.1 7
./wlsh NUS 1.5 50 5 10 0.1 7.5

./wlsh audio 1.5 50 5 10 0.1

./wlsh mnist 1.5 50 5 10 0.1


./wlsh deep1m 1.5 50 5 10 0.1 0.2 #optimal

./wlsh gist 1.5 50 5 10 0.1 0.3 #optimal

./wlsh gist2 1.5 50 5 12 0.18 0.1
./wlsh gist4 1.5 50 5 12 0.16 0.15
./wlsh gist6 1.5 50 5 12 0.14 0.25
./wlsh gist8 1.5 50 5 12 0.12 0.3
./wlsh gist 1.5 50 5 12 0.12 0.3
cd .

./wlsh gist2 1.5 50 5 12 0.15 0.1
./wlsh gist4 1.5 50 5 12 0.13 0.15
./wlsh gist6 1.5 50 5 12 0.12 0.25
./wlsh gist8 1.5 50 5 12 0.11 0.3
./wlsh gist 1.5 50 5 12 0.12 0.3
cd .

./wlsh gist 1.5 50 5 12 -5 0.3 #optimal
./wlsh gist 1.5 50 5 12 -20 0.3 #optimal
cd .

./wlsh tiny 1.5 50 5 12 -5 0.2
./wlsh tiny 1.5 50 5 12 -50 0.2
cd ..

./wlsh sift10M 1.5 50 5 12 0.06 80


./wlsh sift10M 1.5 50 5 12 0.02 30

./wlsh sift100M 1.5 50 5 12 0.04 80

./wlsh Trevi 1.5 50 5 10 -5 600

./wlsh sift10M 1.5 50 5 12 -5 30

./wlsh Trevi 1.5 50 5 12 -5 800 #800 optimal