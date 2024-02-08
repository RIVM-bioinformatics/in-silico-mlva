import argparse, os.path, csv, textwrap, glob
import pandas as pd
import numpy as np
from pathlib import Path
from termcolor import colored
from tqdm import tqdm

def getmylogo(pth):
    exec_globals = {}
    with open(pth, 'r') as lfile:
        exec(lfile.read(), exec_globals)
    logo = exec_globals.get('logo', None)
    return logo

def parse_arguments(logo):
    arg = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent(f"""
        {colored(logo, 'red', attrs=["bold"])}
        {colored('In silico MLVA typing for MRSA:', 'white', attrs=["bold", "underline"])}

        Searches possible products formed based on the primer blast output.
        Script will then determine if these are within a known bin size for each VNTR.
-----------------------------------------------------------------------------------
        {colored('Example usage:', 'green', attrs=["bold", "underline"])}
            python {os.path.abspath(__file__)} 
            --blast_primer example/output_blastn/RIVM_M096462_primers-blastn.csv
            --blast_repeat example/output_blastn/RIVM_M096462_repeat-blastn.csv 
            --output example/output_mlva_typing/)
-----------------------------------------------------------------------------------
        """))
    arg.add_argument(
        "-bp",
        "--blast_primer",
        metavar="Name",
        help="CSV file with primer blastn output",
        type=str,
        required=False,
    )

    arg.add_argument(
        "-br",
        "--blast_repeat",
        metavar="Name",
        help="CSV file with blastn repeat sequence output",
        type=str,
        required=False,
    )

    arg.add_argument(
        "-i",
        "--input",
        metavar="Name",
        help="The original fasta input directory",
        type=str,
        required=True,
    )

    arg.add_argument(
        "-o",
        "--output",
        metavar="Name",
        help="Output directory where you want to copy to",
        type=str,
        required=False,
    )

    return arg.parse_args()

def determine_outdir(flg_out):
    if flg_out == None:
        outdir = os.path.abspath('') + '/mlva_typing'
        Path(os.path.abspath(f"{outdir}")).mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {outdir}")
    else:
        Path(os.path.abspath(f"{flg_out}/mlva_typing")).mkdir(parents=True, exist_ok=True)
        outdir = os.path.abspath(flg_out) + '/mlva_typing'
        print(f"Output directory: {outdir}")
    return outdir

def csv_to_list(file):
    with open(file) as f:
        if os.path.getsize(file) > 0:
            reader = csv.reader(f)
            data = list(reader)
            return data

def determine_repeats_inrange(for1, for2, repeats): # This returns 99 if no primers are found
    repeat_sequence_in_range_fc = []
    for f in for1:
        for seq in repeats:
            if int(f)-int(seq.split('-')[0]) > 0 and int(f)-int(seq.split('-')[0]) <= 1200:
                repeat_sequence_in_range_fc.append(seq)
    for f in for2:
        for seq in repeats:
            if int(f)-int(seq.split('-')[0]) <= 0 and int(f)-int(seq.split('-')[0]) >= -1200:
                repeat_sequence_in_range_fc.append(seq)
    return repeat_sequence_in_range_fc

def determine_chain(repeats):
    sorted_range_list = sorted(repeats, key=lambda x: int(x.split('-')[0])) # Sort all of the ranges with the lowest value to highest starting value
    if len(sorted_range_list) > 0:
        repeat_counter = 1
        for i in range(len(sorted_range_list) - 1):
            current_element = sorted_range_list[i]
            next_element = sorted_range_list[i + 1]
            current_end = int(current_element.split('-')[1])
            next_start = int(next_element.split('-')[0])
            if current_end < next_start: #This is where you keep counting if the end matches the next start.
                if int(current_element.split('-')[1]) + 1 == int(next_element.split('-')[0]):
                    repeat_counter += 1
        return str(repeat_counter)
    else:
        return '99'

def get_number_repeats(dataframe, dataframe_repeat, fname, rname, bitscori, bitscori_r):
    dataframe = dataframe.astype({'bitscore':'float'})
    dataframe = dataframe.loc[dataframe['bitscore'] >= bitscori]
    dataframe_repeat = dataframe_repeat.astype({'bitscore':'float'})
    dataframe_repeat = dataframe_repeat.loc[dataframe_repeat['bitscore'] >= bitscori_r]
    forward_list = list(set(dataframe.loc[dataframe['sseqid'] == fname]['qseqid'].tolist()))
    shared_list = list(set(dataframe.loc[(dataframe['sseqid'] == rname) & (dataframe['qseqid'].isin(forward_list))]['qseqid'].tolist())) # contigs which have results for both forward and reverse
    repeat_list = list(set(dataframe_repeat.loc[dataframe_repeat['sseqid'] == 'VNTR63_01']['qseqid'].tolist()))
    if len(shared_list) == 0: # This likely means that 1 of the primers wasn't found! Which is expected for 63_01 reverse
        if fname != 'VNTR63_01_Ff':
            print(f"Only 1 primer found for {fname}, can't continue")
            exit()
        else:
            if len(forward_list) == 0:
                print(f"Not even forward primer found for {fname}, can't continue")
                exit()
            else:
                forward_pos1, forward_pos2 = [], []
                for l in forward_list:
                    forward_pos1_pre = dataframe.loc[(dataframe['qseqid'] == l) & (dataframe['sseqid'] == fname)]['qend'].tolist()
                    forward_pos2_pre = dataframe.loc[(dataframe['qseqid'] == l) & (dataframe['sseqid'] == fname)]['qstart'].tolist()
                    forward_pos1.extend(forward_pos1_pre)
                    forward_pos2.extend(forward_pos2_pre)
    else:
        forward_pos1, forward_pos2 = [], []
        for l in shared_list: # for every contig
            forward_pos1_pre = dataframe.loc[(dataframe['qseqid'] == l) & (dataframe['sseqid'] == fname)]['qend'].tolist()
            forward_pos1.extend(forward_pos1_pre)
            forward_pos2_pre = dataframe.loc[(dataframe['qseqid'] == l) & (dataframe['sseqid'] == fname)]['qstart'].tolist()
            forward_pos2.extend(forward_pos2_pre)
    repeat_pos_list = []

    for l in repeat_list:
        repeat_pos = dataframe_repeat.loc[(dataframe_repeat['qseqid'] == l) & (dataframe_repeat['sseqid'] == 'VNTR63_01')].apply(lambda row: f"{row['qstart']}-{row['qend']}", axis=1).tolist()
        repeat_pos_list.extend(repeat_pos)
    repeat_sequence_in_range = determine_repeats_inrange(forward_pos1, forward_pos2, repeat_pos_list) # Only these make sense because they are in range with the forward primer
    number_of_repeats = determine_chain(repeat_sequence_in_range)
    return number_of_repeats

def get_possible_sizes(dataframe, fname, rname, bitscori):
    # TODO: Remove any possible sizes that are not within their found repeat locations.
    forward_list = list(set(dataframe.loc[dataframe['sseqid'] == fname]['qseqid'].tolist()))
    shared_list = list(set(dataframe.loc[(dataframe['sseqid'] == rname) & (dataframe['qseqid'].isin(forward_list))]['qseqid'].tolist())) # contigs which have results for both forward and reverse
    all_posi_fc = []
    dataframe = dataframe.astype({'bitscore':'float'})
    dataframe = dataframe.loc[dataframe['bitscore'] >= bitscori]
    for l in shared_list: # for every contig
        forward_pos1 = dataframe.loc[(dataframe['qseqid'] == l) & (dataframe['sseqid'] == fname)]['qend'].tolist()
        reverse_pos1 = dataframe.loc[(dataframe['qseqid'] == l) & (dataframe['sseqid'] == rname)]['qstart'].tolist()
        forward_pos2 = dataframe.loc[(dataframe['qseqid'] == l) & (dataframe['sseqid'] == fname)]['qstart'].tolist()
        reverse_pos2 = dataframe.loc[(dataframe['qseqid'] == l) & (dataframe['sseqid'] == rname)]['qend'].tolist()
        for f in forward_pos1: # Because I don't know the orientation of the chromosome the forward might actually be at a higher location than the reverse. So by subtracting reverse from forward you should have a positive integer <1200 whenever this happens.
            for r in reverse_pos1: # Taking the qend of forward en qstart of reverse means you calculate the total size of the product whenever the forward is located upstream compared to the reverse.
                if int(f)-int(r) >= 0 and int(f)-int(r) <= 1200:
                    all_posi_fc.append(int(f)-int(r))
        for f in forward_pos2:
            for r in reverse_pos2:
                if int(f)-int(r) >= -1200 and int(f)-int(r) <= 0:
                    all_posi_fc.append(0-(int(f)-int(r)))
    return list(set(all_posi_fc))

def get_mlva_dict(dataframe):
    MLVA_dict_fc = {}
    MLVA_dict_fc['MLVA_MecA'] = get_possible_sizes(dataframe,'MLVA_MecA_Ff','MLVA_MecA_r',30)
    MLVA_dict_fc['VNTR09_01'] = get_possible_sizes(dataframe,'VNTR09_01_Ff','VNTR09_01_r',30)
    MLVA_dict_fc['VNTR61_01'] = get_possible_sizes(dataframe,'VNTR61_01_Nf','VNTR61_01_r',30)
    MLVA_dict_fc['VNTR61_02'] = get_possible_sizes(dataframe,'VNTR61_02_Vf','VNTR61_02_r',30)
    MLVA_dict_fc['VNTR67_01'] = get_possible_sizes(dataframe,'VNTR67_01_Pf','VNTR67_01_r',30)
    MLVA_dict_fc['MLVA_PVL'] = get_possible_sizes(dataframe,'MLVA_PVL_Ff','MLVA_PVL_r',30)
    MLVA_dict_fc['VNTR21_01'] = get_possible_sizes(dataframe,'VNTR21_01_Vf','VNTR21_01_r',30)
    MLVA_dict_fc['VNTR24_01'] = get_possible_sizes(dataframe,'VNTR24_01_Pf','VNTR24_01_r',30)
    MLVA_dict_fc['VNTR63_01'] = get_possible_sizes(dataframe,'VNTR63_01_Ff','VNTR63_01_r',15)
    MLVA_dict_fc['VNTR81_01'] = get_possible_sizes(dataframe,'VNTR81_01_Nf','VNTR81_01_r',30)
    return MLVA_dict_fc

def mec_or_pvl(df_mappings, mecpvl_list, mlvadict):
    mp_dict_fc = {}
    for mp in mecpvl_list:
        values = mlvadict[mp]
        if len(values) == 0:
            if mp == 'MLVA_MecA':
                mp_dict_fc[mp] = "MecA MecC Negative"
            elif mp == 'MLVA_PVL':
                mp_dict_fc[mp] = "PVL Negative"
        elif len(set(values)) == 1:
            try:
                value = df_mappings.loc[(df_mappings['VNTR'] == mp) & (df_mappings['Start'] <= values[0]) & (df_mappings['Stop'] >= values[0])].Value.item()
                if mp == 'MLVA_MecA':
                    if value == 1:
                        mp_dict_fc[mp] = "MecA Positive"
                    elif value == 2:
                        mp_dict_fc[mp] = "MecC Positive"
                elif mp == 'MLVA_PVL':
                    if value == 1:
                        mp_dict_fc[mp] = "PVL Positive"
            except:
                mp_dict_fc[mp] = "Something failed"
    return mp_dict_fc

def get_my_profile(df_mappings, vntr_lst, mlvadict, df, df2):
    profile_fc = "" 
    profile_fc_list = []
    deviated = False
    for v in vntr_lst:
        profile_fc = profile_fc.strip('-') # should just not add a '-' on the first value but this works too
        if v == 'VNTR63_01': # Making an exception for VNTR63_01 because the reverse primer can't be found in about 30% of the isolates.
            repeat_found = get_number_repeats(df,df2,'VNTR63_01_Ff','VNTR63_01_r',25,55)
            profile_fc = f"{profile_fc}-{repeat_found}"
        else:
            values = mlvadict[v] # These values should be checked if they are within range of their primers!
            if deviated == False:
                if len(values) == 0:
                    profile_fc = f"{profile_fc}-99"
                elif len(values) == 1:
                    try: # vntr_no = df_mappings.loc[(df_mappings['VNTR'] == v) & (df_mappings['Start'] <= values[0]) & (df_mappings['Stop'] >= values[0])].Value.item()
                        matching_ranges = df_mappings.loc[
                        (df_mappings['VNTR'] == v) & 
                        (df_mappings['Start'] <= values[0]) & 
                        (df_mappings['Stop'] >= values[0])]
                        if not matching_ranges.empty:
                            vntr_no = matching_ranges['Value'].item()
                        else: # If the found size is outside a range it takes the closest value, perhaps should print a message for an isolate when this has happened.
                            closest_value_index = np.argmin(
                            np.abs(df_mappings.loc[df_mappings['VNTR'] == v, 'Start'].values - values[0]) +
                            np.abs(df_mappings.loc[df_mappings['VNTR'] == v, 'Stop'].values - values[0]))
                            vntr_no = df_mappings.loc[df_mappings['VNTR'] == v, 'Value'].iloc[closest_value_index]
                        profile_fc = f"{profile_fc}-{vntr_no}"
                    except:
                        profile_fc = f"{profile_fc}-99"
                else: # This means multiple locations have been found, however they might be in the same bin so this checks if that's the case or not.
                    bin_list = []
                    for index, val in enumerate(values):
                        try:
                            vntr_no = df_mappings.loc[(df_mappings['VNTR'] == v) & (df_mappings['Start'] <= val) & (df_mappings['Stop'] >= val)].Value.item()
                            bin_list.append(vntr_no)
                        except:
                            bin_list.append('99')
                        if index == len(values) - 1:
                            if len(list(set(bin_list))) != 1:
                                deviated = True
                            else:
                                profile_fc = f"{profile_fc}-{bin_list[0]}"
                    if deviated == True:
                        possibles = []
                        for n in values:
                            try:
                                vntr_no = df_mappings.loc[(df_mappings['VNTR'] == v) & (df_mappings['Start'] <= n) & (df_mappings['Stop'] >= n)].Value.item()
                                possibles.append(vntr_no)
                            except:
                                possibles.append('99')
                        for p in possibles:
                            profile_fc_list.append(f"{profile_fc}-{p}")
            else: # This happens if deviates = True and multiple profiles were found, can probably be excluded because 
                deviated_fc = determine_deviated_profiles(profile_fc_list, values, df_mappings, v)
    if len(profile_fc_list) == 0: # This is done to not add a profile when multiples were found but can probably remove because it shouldn't happen. Besides how would we determine which one would be correct
        profile_fc_list.append(profile_fc)
    profile_fc_list = double_pad(list(set(profile_fc_list)))
    return profile_fc_list

def determine_deviated_profiles(a_profile_list, val, dfmap, fcv):
    times_to_pop = len(a_profile_list)
    for pli in range(0,len(a_profile_list)):
        pl = a_profile_list[pli]
        if len(val) == 0:
            pl = f"{pl}-99"
            a_profile_list.append(pl)
        elif len(val) == 1:
            try:
                vntr_no = dfmap.loc[(dfmap['VNTR'] == fcv) & (dfmap['Start'] <= val[0]) & (dfmap['Stop'] >= val[0])].Value.item()
                pl = f"{pl}-{vntr_no}"
                successful = True
            except:
                pl = f"{pl}-99"
                successful = False
                a_profile_list.append(pl)
            if successful:
                a_profile_list.append(pl)
        else:
            possibles = []
            for n in val:
                try:
                    vntr_no = dfmap.loc[(dfmap['VNTR'] == fcv) & (dfmap['Start'] <= n) & (dfmap['Stop'] >= n)].Value.item()
                    possibles.append(vntr_no)
                except:
                    possibles.append('99')
            for p in possibles:
                a_profile_list.append(f"{pl}-{p}")
    for p in range(0,times_to_pop):
        a_profile_list.pop(0)
    return a_profile_list

def double_pad(profile_lst):
    double_padded_lst_fc = []
    for l in profile_lst:
        ls = l.split('-')
        new_l_fc = f""
        for n in ls:
            if len(n) == 1:
                new_l_fc = f"{new_l_fc}-0{n}".strip('-')
            else:
                new_l_fc = f"{new_l_fc}-{n}".strip('-')
        double_padded_lst_fc.append(new_l_fc)
    return double_padded_lst_fc

def write_to_file(profile_lst, df_mappings, mecpvl_list, mlvadict, fname, outd):
    base_input_fc = os.path.basename(fname[0]).replace('_primers-blastn.csv','')
    with open(f"{outd}/{base_input_fc}_MLVA.txt", "w") as my_file:
        # my_file.write(str(profile_lst) + '\n')
        for p in profile_lst:
            output_mecpvl = mec_or_pvl(df_mappings,mecpvl_list,mlvadict)
            my_file.write(f"MLVA profile: {p}" + '\n')
            my_file.write(output_mecpvl['MLVA_MecA'] + '\n')
            my_file.write(output_mecpvl['MLVA_PVL'] + '\n')

def main():
    current_file_path = os.path.abspath(__file__)
    parent_dir_path = os.path.dirname(os.path.dirname(current_file_path))
    mrsa_mapping = os.path.join(parent_dir_path, "files", "mrsa_mappings.csv")
    logo_path = os.path.join(parent_dir_path, "files", "logo.txt")
    df_mapping = pd.read_csv(mrsa_mapping, sep=",")

    flags = parse_arguments(getmylogo(logo_path))
    outdir = determine_outdir(flags.output)

    blastn_header = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']
    static_list = ['MLVA_MecA', 'MLVA_PVL']
    vntr_list = ['VNTR09_01', 'VNTR61_01', 'VNTR61_02', 'VNTR67_01', 'VNTR21_01', 'VNTR24_01', 'VNTR63_01', 'VNTR81_01']
    list_of_files = glob.glob(os.path.abspath(f"{flags.input}/*"))
    for file in list_of_files:
        basename = os.path.splitext(os.path.basename(file))[0]
        outputname = f"{outdir}/blastn/{basename}"
        # blast_input = [flags.blast_primer]
        blast_input = [f"{outputname}_primers-blastn.csv"]
        # repeat_file = [flags.blast_repeat]
        repeat_file = [f"{outputname}_repeat-blastn.csv"]
        single_entry_list_primers = [entry for file in (csv_to_list(f) for f in blast_input) if file for entry in file]
        single_entry_list_repeats = [entry for file in (csv_to_list(f) for f in repeat_file) if file for entry in file]
        df = pd.DataFrame(single_entry_list_primers, columns=blastn_header) # The blast primer output to a df
        df2 = pd.DataFrame(single_entry_list_repeats, columns=blastn_header) # The blast repeat output to a df

        MLVA_dict = get_mlva_dict(df)
        profiles_in_a_list = get_my_profile(df_mapping, vntr_list, MLVA_dict, df, df2)

        in_silico_profile = profiles_in_a_list[0]

        write_to_file(profiles_in_a_list,df_mapping,static_list,MLVA_dict,blast_input,outdir)
        print(f"profile for {file}: {in_silico_profile}")

if __name__ == "__main__":
    main()
