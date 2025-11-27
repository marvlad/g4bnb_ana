source /al9g4bnb/local/root_install/bin/thisroot.sh
source /al9g4bnb/local/geant4-v10.6.1_install/bin/geant4.sh

export BOOST_LIB=/al9g4bnb/local/boost_1_85_0_install/lib
export BOOST_INCLUDE=/al9g4bnb/local/boost_1_85_0_install/include

export CPLUS_INCLUDE_PATH=/al9g4bnb/local/boost_1_85_0_install/include:${CPLUS_INCLUDE_PATH}
export LD_LIBRARY_PATH=/al9g4bnb/local/boost_1_85_0_install/lib:$LD_LIBRARY_PATH

#export DK2NU=/al9g4bnb/local/dk2nu_install
export DK2NU=/exp/annie/app/users/mascenci/g4bnb_detached/dk2nu/install
echo $DK2NU
export DK2NUDATA_INC=$DK2NU/include
export DK2NUDATA_LIB=$DK2NU/lib
export LD_LIBRARY_PATH=$DK2NUDATA_LIB:$LD_LIBRARY_PATH
