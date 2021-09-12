# Consistent LCG environment (http://lcginfo.cern.ch)
release=`cat /etc/redhat-release`
if test "${release#*release 6.}" != "$release"; then
  # Scientific Linux 6
  . /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_97python3 x86_64-slc6-gcc8-opt
  export SCRAM_ARCH=slc6_amd64_gcc830
else
  # CentOS 7
  . /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_97python3 x86_64-centos7-gcc9-opt
  export SCRAM_ARCH=slc7_amd64_gcc920
fi

export CMAKE_BINARY_DIR=$(pwd)
export PATH=$CMAKE_BINARY_DIR/bin:$PATH
