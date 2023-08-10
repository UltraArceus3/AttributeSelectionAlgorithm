#!/bin/bash

set -e

directory_json="/home/nachiket/RLA_CL_EXTRACT"

cd -- "$directory_json"

make

bin/rlacl in/config/config_q_gram_ds2.xml  100000 0 0 0 0
bin/rlacl in/config/config_q_gram_ds3.xml  200000 0 0 0 0
bin/rlacl in/config/config_q_gram_ds4.xml  400000 0 0 0 0
bin/rlacl in/config/config_q_gram_ds5.xml  600000 0 0 0 0
bin/rlacl in/config/config_q_gram_ds6.xml  800000 0 0 0 0


# bin/rlacl in/config/config_edit_ds2.xml  100000 0 0 0 0
# bin/rlacl in/config/config_edit_ds3.xml  200000 0 0 0 0
# bin/rlacl in/config/config_edit_ds4.xml  400000 0 0 0 0
# bin/rlacl in/config/config_edit_ds5.xml  600000 0 0 0 0
# bin/rlacl in/config/config_edit_ds6.xml  800000 0 0 0 0







# bin/rlacl in/config.xml 100000 0 0 0 0
# bin/rlacl in/config.xml 200000 0 0 0 0
# bin/rlacl in/config.xml 300000 0 0 0 0
# bin/rlacl in/config.xml 400000 0 0 0 0
# bin/rlacl in/config.xml 500000 0 0 0 0
# bin/rlacl in/config.xml 600000 0 0 0 0
# bin/rlacl in/config.xml 700000 0 0 0 0
# bin/rlacl in/config.xml 800000 0 0 0 0
# bin/rlacl in/config.xml 900000 0 0 0 0
# bin/rlacl in/config.xml 1000000 0 0 0 0

# bin/rlacl in/config.xml 100000 701 999983 99991 900001
# bin/rlacl in/config.xml 200000 701 999983 99991 900001
# bin/rlacl in/config.xml 300000 701 999983 99991 900001
# bin/rlacl in/config.xml 400000 701 999983 99991 900001
# bin/rlacl in/config.xml 500000 701 999983 99991 900001
# bin/rlacl in/config.xml 600000 701 999983 99991 900001
# bin/rlacl in/config.xml 700000 701 999983 99991 900001
# bin/rlacl in/config.xml 800000 701 999983 99991 900001
# bin/rlacl in/config.xml 900000 701 999983 99991 900001
# bin/rlacl in/config.xml 1000000 701 999983 99991 900001

# bin/rlacl in/config.xml 100000 6607 99991 58211 55001
# bin/rlacl in/config.xml 200000 6607 99991 58211 55001
# bin/rlacl in/config.xml 300000 6607 99991 58211 55001
# bin/rlacl in/config.xml 400000 6607 99991 58211 55001
# bin/rlacl in/config.xml 500000 6607 99991 58211 55001
# bin/rlacl in/config.xml 600000 6607 99991 58211 55001
# bin/rlacl in/config.xml 700000 6607 99991 58211 55001
# bin/rlacl in/config.xml 800000 6607 99991 58211 55001
# bin/rlacl in/config.xml 900000 6607 99991 58211 55001
# bin/rlacl in/config.xml 1000000 6607 99991 58211 55001


