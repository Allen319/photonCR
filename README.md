# photonCR-analyzer
This repository uses Git and hence should be cloned with
```sh
git clone https://github.com/Allen319/photonCR.git
```
No CMSSW environment is needed, so it is unnecessary to put this repository under CMSSW.

We use [`LCG`](http://lcginfo.cern.ch)  on lxplus.
At the start of each session, set up the environment with

```sh
. ./LCG_env.sh
```
If you are running locally, you will need libraries as following
`yaml-cpp`, `boost` and `ROOT 6`.
These environments should be already setup on your local PC. 
And then set up the environment with

```sh
. ./local_env.sh
```
Build the package with the following commands:

```sh
mkdir build
cd build
cmake ..
make -j $(nproc)
cd -
```

Simply run with 
```sh
mkdir -p outputs/DY
mkdir outputs/GJet
skim -c closure_new.yaml
runphoton -c step2.yaml
```
