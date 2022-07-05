import subprocess
import sys
import glob
import os
import re
import gzip
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from datetime import datetime
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import factorial
import matplotlib
from numpy import trapz
from math import gamma

matplotlib.style.use('ggplot')
np.seterr(over='ignore')
flag_for_asm = 0
# measure time
startTime = datetime.now()
print(str(datetime.now() - startTime) + " starting")

n_runs = 0
past_runs = int(sys.argv[3]) + 1
number_of_input_file = sys.argv[4]
path_to_lists = sys.argv[5]


def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1
        yield tuple(pool[i] for i in indices)


def get_fasta_of_n_random_genomes():  # Require /mnt/g/from_dust/genbank/bacteria
    number_of_genomes_ready = 0
    global glob_target_genomes
    # glob_target_genomes = random.sample(genomes_list,n) # was in first place
    with open(path_to_lists + "/a_list_" +
              number_of_input_file + ".txt", "r") as lets_do_them_all:
        glob_target_genomes = lets_do_them_all.readline().split(",")[0].split("_and_")
    for refseq in glob.glob('/mnt/g/from_dust/genbank/bacteria/**/*',
                            recursive=True):  # dir can be moved
        if refseq[-3:] == '.gz' and refseq[-7:] != 'D5SUMS/' and len(glob_target_genomes) >= 2:
            for k, sampl in enumerate(glob_target_genomes):
                print(refseq.split('_'))
                if refseq.split('_')[3] == sampl:  # fix this _ shit, asshole
                    number_of_genomes_ready += 1
                    with gzip.open(refseq, 'rb') as inf:  # gzip in bytes
                        e = inf.read()
                    with open(path + 'genomes/trash.txt', 'w') as f2:
                        f2.write('\n' + e.decode('utf-8'))
                    with open(path + 'genomes/trash.txt', 'r') as f3:
                        s = []
                        t = 0
                        for line in f3:
                            if line[0] == '>':
                                t += 1  # counter, cause we only compute complete genomes
                                if t == 2:
                                    break
                            if ">" in line:
                                s.append(">" + line.split(">")[1])
                            else:
                                s.append(line)
                    with open(path + 'genomes/{0}_{1}'.format(sampl, k), 'a') as f1:
                        f1.write(''.join(s))
                    with open(path + 'genomes/{0}_{1}'.format(sampl, k), 'r') as inp:
                        with open(path + 'genomes/{0}_{1}.fasta'.format(sampl, k), 'w') as outp:
                            for line_number, line in enumerate(inp):
                                if line_number != 0:
                                    outp.write(line)
                    os.remove(path + 'genomes/{0}_{1}'.format(sampl, k))
    try:
        os.remove(path + 'genomes/trash.txt')
    except FileNotFoundError:
        pass


def mauve_files(directory, out):  # dir with fastas in / mauve out file
    right_call = ['progressiveMauve', '--output={0}blocks/my_seqs.xmfa'.format(path)]
    for filename in glob.glob(directory + "/*"):
        right_call.append(filename)
    print(right_call)
    subprocess.call(right_call)
    os.system("cp {0}blocks/my_seqs.xmfa {1}fast_operations".format(path, path))


def mauve_to_blocks(filename, my_dir):  # filename is mauve file / my_dir is where
    p = re.compile(u"=\n")
    prod = re.compile(u">.*\n")
    with open(my_dir + filename, "r") as f:
        result = ""
        current = 0
        for line in f:
            result += line
            if re.match(prod, line):
                current += 1
            if re.match(p, line):
                with open(my_dir + filename + "__" + str(current), "a") as f2:
                    f2.write(result)
                    f2.close()
                current = 0
                result = ""

    usiversal_max = 0
    for filename in glob.glob(my_dir + "*", recursive=True):
        if "__" in filename:
            if int(filename.split("__")[-1]) > usiversal_max:
                usiversal_max = int(filename.split("__")[-1])
    for every_file in glob.glob(my_dir + "*", recursive=True):
        if str(every_file[-len(str(usiversal_max)):]) != str(usiversal_max):
            try:
                os.remove(every_file)
            except IsADirectoryError:
                pass


def big_mauve_blocks_to_thousands(directory):
    p = re.compile(u"=")  # end of big block
    prod = re.compile(u">")  # begining of small block
    A = re.compile(u"A")
    T = re.compile(u"T")
    C = re.compile(u"C")
    G = re.compile(u"G")
    gap = re.compile(u"-")
    abz = re.compile(u"\n")
    for filename in glob.glob(directory + "*"):
        with open(filename, "r") as f:
            block_number = 0
            file_counter = 0
            align = 0
            buff = ""
            for line in f:
                if re.match(p, line):
                    block_number += 1
                    file_counter = 0
                    buff = ""
                    align = 0
                elif re.match(prod, line):
                    file_counter = 0
                    buff = ""
                    align += 1
                elif (re.match(A, line) or re.match(C, line) or re.match(G, line) or re.match(T, line) or
                      re.match(gap, line)):
                    for c in line:
                        buff += c
                        if (len(buff) == 1012 and len(buff.split('\n')) == 13) or (
                                len(buff) == 1013 and len(buff.split('\n')) == 14) or (
                                len(buff) == 1014 and len(buff.split('\n')) == 15) or (
                                len(buff) == 1015 and len(buff.split('\n')) == 16):
                            with open(directory + str(block_number) + "_" + str(file_counter), "a") as f2:
                                f2.write("\n> " + str(align) + "_" + str(block_number) + "_" + str(file_counter)
                                         + "\n" + buff + "\n=")
                                buff = ""
                                file_counter += 1


def third_script(directory):  # clean block directory from '-'
    for filename in glob.glob(directory + "/**/*", recursive=True):
        try:
            with open(filename, "r") as f:
                count = 0
                for line in f:
                    match1 = re.findall(u"--", line)
                    count += len(match1)
            if count >= 1:
                os.remove(filename)
        except IsADirectoryError:
            continue


def delete_all_files_in_dir(directory):
    for every_file in glob.glob(directory + "/**"):
        try:
            os.remove(every_file)
        except IsADirectoryError:
            continue


def fourth_script(directory):
    asms = []  # block for number of genomes
    my_counter = 0
    for filename in glob.glob(path + 'blocks' + "/**/*", recursive=True):
        if my_counter == 0:
            my_counter += 1
            with open(filename, "r") as f:
                for line in f:
                    if line[0] == ">":
                        asms.append(line)

    n_genomes = len(asms)
    realy_list = []
    sum_n = 0
    counter_gen_numb = n_genomes
    while counter_gen_numb > 0:
        sum_n += counter_gen_numb
        counter_gen_numb = counter_gen_numb - 1

    # creating list for blocks
    new_list = [[i] for i in range(0, sum_n)]
    for item in new_list:
        item.pop(0)

    for filename in glob.glob(directory + 'blocks' + "/**/*", recursive=True):
        c = 0
        f = open(filename, "r", encoding="utf-8")
        x = f.readlines()
        my_str = ''.join(x)
        final = my_str.split("=")
        if "\n" in final:
            final.remove("\n")
        if "" in final:
            final.remove("")
        f.close()
        final = tuple(final)
        for first in final:
            for second in final:
                if final.index(second) >= final.index(first):
                    if first == second:
                        new_list[c].append(0)
                        c += 1
                    else:
                        with open(path + 'gr' + "/Doc0", "w", encoding="utf-8") as f2:  # can be optimized here
                            f2.write(first + second)
                        with open(path + 'gr' + "/Doc0", "r") as inp:
                            with open(path + 'gr' + "/Doc1", "w") as outp:
                                for k, line in enumerate(inp):
                                    if k != 0:
                                        outp.write(line)
                        subprocess.call(['snp-sites', '-p', '-o', path + 'gr' + "/Doc2", path + 'gr' +
                                                   "/Doc1"])
                        try:
                            with open(path + 'gr' + "/Doc2", "r", encoding="utf-8") as f3:
                                for numb_of_lin, lines in enumerate(f3):
                                    if numb_of_lin == 0:
                                        new_list[c].append(int(lines.split(" ")[1][:-1]))
                        except FileNotFoundError:
                            new_list[c].append(0)
                        try:
                            os.remove(path + 'gr' + "/Doc2")
                        except FileNotFoundError:
                            pass
                        os.remove(path + 'gr' + "/Doc0")
                        os.remove(path + 'gr' + "/Doc1")
                        c += 1
    another_list = [Counter(thing) for thing in new_list]

    e = [i for i in range(1, n_genomes + 1)]

    for my_numb in e:
        while realy_list.count(my_numb) != my_numb:
            realy_list.append(my_numb)
    realy_list.reverse()

    # for many empty lists
    many_lists = [i for i in range(1, 1000)]
    for i in many_lists:
        exec("lis{0} = []".format(str(i)))

    for k, j in enumerate(realy_list):  # it was huge construction
        exec("lis{0}.append(another_list[k])".format(str(n_genomes + 1 - j)))
    """
    # u can prewrite "lis" from previous runs here and comment above to skip this steps (remember about counter below)
    """
    for i in range(1, n_genomes + 1):
        with open(path + "lists", "a") as f4:
            exec("f4.write(\"lis{0}=\".format(i))".format(i))
            exec("f4.write(\"[\"+ \",\".join(" + "str(v) for v in " + "lis{0}".format(i) + ")+\"]\"+\"\\n\")")

    # remember metadata
    with open(path + "names", "a") as f5:
        f5.write(glob_target_genomes[0] + "\n")
        f5.write(glob_target_genomes[1] + "\n")
    with open(path + "start_pars", "a") as f6:
        f6.write(sys.argv[1] + "\n" + sys.argv[2] + "\n" + sys.argv[3] + "\n")

    for group in combinations(range(1, n_genomes + 1), 2):
        f_t = [[] for i in range(1000)]
        exec("exec1 = lis{0}[{1}].items()".format(str(group[0]), str(group[1] - group[0])))
        for i in exec1:
            f_t[i[0]].append(i[1])
        for k, el in enumerate(f_t):
            f_t[k] = sorted(el)

        fxd = []
        rdyList = {}
        print(f_t)
        for lisT in f_t:
            if lisT == []:
                fxd.append([0])
            else:
                fxd.append(lisT)
        first_zero = fxd.index([0])
        with open(path + "/fast_operations/first_zero.txt", "w") as first_zero_file:
            first_zero_file.write(str(first_zero) + "\n")
            fxd_to_write = []
            for el in fxd:
                fxd_to_write.append(str(el[0]))
            first_zero_file.write(",".join(fxd_to_write) + "\n")
        while len(fxd) >= first_zero:
            fxd.pop()
            if len(fxd) < 2:
                break

        # warning about small list (stains are too close)
        if len(fxd) < 2:
            print("problems in " + str(group[0]) + " " + str(group[1]))
            print(fxd)
            fxd.append([1])
            fxd.append([1])
        for k, el in enumerate(fxd):
            average = sum(el) / len(el)
            rdyList[k] = average
        thing = pd.Series(rdyList)
        with open(path + 'gr' + "/" + str(group[0]) + " " + str(group[1]) + ".txt", "w") as f:
            f.write(str(rdyList))

        plt.xlabel("Number of SNPs")
        plt.ylabel("Blocks with SNPs")
        thing.plot(x="thing.index", y="thing.values")
        plt.ylim(0, 200)
        plt.xlim(0, 100)
        plt.savefig(path + 'gr' + "/" + "single double " + str(group[0]) + " " + str(group[1]) + '.png', format='png',
                    dpi=100)
        plt.clf()


def func(x, a, l, z):
    return ((a ** x) * np.exp(-a) / factorial(x)) * z + (x ** (k - 1) * np.exp(-x / l) * (1 - z)) / (
            (l ** k) * gamma(k))


def fifth_script(path):
    pars = [x * 0.2 for x in range(6, 50)]

    for filename in glob.glob(os.path.join(path + 'gr' + "/", '*.txt')):
        name = filename.split("/")[-1]
        with open(filename, "r") as f:  # might be some problems with "counter" in file
            mf = f.read()
            mf = eval(mf)
            c = 0
            ef = {}
            for val in mf:
                c += mf[val]
            for val in mf:
                ef[val] = mf[val] / c
            ready_lists = []
            for a in ef.values():
                ready_lists.append(a)
            # Compute the area using the composite trapezoidal rule.
            area = trapz(ready_lists, dx=1)
            true_area = area * 0.2
            whole_list = ready_lists.copy()
            while area > true_area:
                ready_lists.pop()
                area = trapz(ready_lists, dx=1)
            minX = len(ready_lists)
            while len(ready_lists) + 1 != len(whole_list):
                if whole_list != []:
                    whole_list.pop()
                else:
                    break
            area_more = trapz(whole_list, dx=1)
            deltaX = (true_area - area) / (area_more - area)
            deltaX = round(deltaX, 2)
            true_X = deltaX + minX
            df = pd.Series(ef) \
                .to_frame('y')
            remember = 1100000000
            for k in pars:
                # def func(x, a, l, z):

                popt, pcov = curve_fit(func, df.index, df['y'], bounds=([0.01, 0.01, 0.01], [5.0, 10.0, 1.0]),
                                       p0=[1, 5, 1.0], max_nfev=10000000)

                # (f(xdata, *popt) - ydata)**2
                sq = sum(((func(df.index, popt[0], popt[1], popt[2]) - df['y']) ** 2) / len(df) * 1000)
                if sq < remember:
                    remember = sq
                    with open(path + 'pars' + "/" + name, "w", encoding="utf-8") as f5:
                        f5.write("name ={0} sq={1} z={2} k={3} a={4} l={5}".format(name, sq, popt[2], k,
                                                                                   popt[0], popt[1]))
                    plt.ylim(0, 0.3)
                    plt.xlim(0, 50)
                    plt.scatter(df.index, df['y'], color='orange')  # scatter or loglog
                    plt.scatter(df.index, func(df.index, *popt))
                    y1, x1 = [0.000001, 0.1], [true_X, true_X]
                    plt.plot(x1, y1, marker='o')
                    plt.title("X = " + str(true_X))
                    plt.savefig(path + 'pars' + "/pars.png", format='png', dpi=100)
                    plt.draw()
                    plt.close("all")


while past_runs < 1000000000:  # look at top
    choosen_philo = sys.argv[2]  # a,b1_1,b2,d1,d2 or e
    # get first pair from list file
    with open(path_to_lists + "/a_list_" +
              number_of_input_file + ".txt", "r") as lets_do_them_all:
        names_of_current_pair = lets_do_them_all.readline().split(",")[0]
    path = sys.argv[1] + '/{0}_{1}_{2}/'.format(str(past_runs), choosen_philo, names_of_current_pair)

    # make some dirs for different steps
    if not os.path.isdir(path + 'blocks'):
        os.makedirs(path + 'blocks')
    if not os.path.isdir(path + 'gr'):
        os.makedirs(path + 'gr')
    if not os.path.isdir(path + 'pars'):
        os.makedirs(path + 'pars')
    if not os.path.isdir(path + 'genomes'):
        os.makedirs(path + 'genomes')
    if not os.path.isdir(path + 'fast_operations'):
        os.makedirs(path + 'fast_operations')

    get_fasta_of_n_random_genomes()
    mauve_files(path + 'genomes', '{0}blocks/my_seqs.xmfa'.format(path))
    mauve_to_blocks('my_seqs.xmfa', path + 'blocks/')
    big_mauve_blocks_to_thousands(path + 'blocks/')
    third_script(path + 'blocks')
    with open(path + "/fast_operations/how_many_blocks.txt", "w") as how_many_blocks:
        how_many_blocks.write(str(len(glob.glob(path + "/blocks/*"))))
    delete_all_files_in_dir(path + 'genomes')
    fourth_script(path + 'genomes')
    fifth_script(sys.argv[1] + '/{0}_{1}_{2}/'.format(str(past_runs), choosen_philo, names_of_current_pair))
    delete_all_files_in_dir(path + 'blocks')

    past_runs += 1
    n_runs += 1
    with open(path_to_lists + "/a_list_" +
              number_of_input_file + ".txt", "r") as old:
        with open(path_to_lists + "/buff_a_list_" +
                  number_of_input_file + ".txt", "w") as new:
            for k, line in enumerate(old):
                if k != 0:
                    new.write(line)
    with open(path_to_lists + "/buff_a_list_" +
              number_of_input_file + ".txt", "r") as old:
        with open(path_to_lists + "/a_list_" +
                  number_of_input_file + ".txt", "w") as new:
            for line in old:
                new.write(line)
    print(str(datetime.now() - startTime) + " done " + str(n_runs))
