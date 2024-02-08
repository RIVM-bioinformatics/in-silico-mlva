# Define some colors and text formatting
BOLD=$(tput bold)
RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
WHITE=$(tput setaf 7)
RESET=$(tput sgr0)

echo " "
echo "${RED}         _             _  _  _            __  __  _  __   __ _   ${RESET}"
echo "${RED}        (_) _ _    ___(_)| |(_) __  ___  |  \/  || | \ \ / //_\  ${RESET}"
echo "${RED}        | || ' \  (_-<| || || |/ _|/ _ \ | |\/| || |__\ V // _ \ ${RESET}"
echo "${RED}        |_||_||_| /__/|_||_||_|\__|\___/ |_|  |_||____|\_//_/ \_\ ${RESET}"
echo "${RED}        ${RESET}"

echo "${BOLD}${WHITE}       In silico MLVA typing for MRSA${RESET}"
echo " "
echo "       Searches possible products formed based on the primer blast output."
echo "       Script will then determine if these are within a known bin size for each VNTR."
echo " "