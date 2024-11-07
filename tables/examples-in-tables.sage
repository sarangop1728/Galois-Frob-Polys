# Return the Weil polynomial corresponding to an LMFDB label.
def weil_poly_from_label(Pol, s):
    l = s.split(".")
    g = int(l[0])
    q = int(l[1])
    l2 = l[2].split("_")
    l3 = [1]
    for s1 in l2:
        sg = -1 if s1[0] == 'a' else 1
        t = [ord(j)-97 for j in s1]
        t.reverse()
        l3.append(sg*ZZ(t, base=26))
    l4 = l3 + [q^i*l3[-i-1] for i in range(1,len(l3))]
    l4.reverse()
    return Pol(l4)

dim2_A_s = ["2.2.ac_d","2.2.ad_f","2.3.ad_f"]
dim2_B_s = ["2.2.ab_a"]
dim2_C_s = ["2.2.ac_c","2.4.ac_e","2.2.a_ae"]
dim2_A_ns = ["2.3.ad_i","2.2.a_d"]
dim2_B_ns = ["2.2.ad_g", "2.4.ah_u"]
dim2_C_ns = ["2.2.ac_e", "2.2.a_a", "2.4.ag_q", "2.4.a_ai"]
dim3_A = ["3.2.ad_f_ah", "3.2.ad_g_aj", "3.2.a_a_ad", "3.2.ac_a_d", "3.2.ae_j_ap", "3.7.ak_bw_afv"]
dim3_B = ["3.2.ab_ab_c", "3.4.ac_ab_g", "3.2.ac_b_a"]
dim3_C = ["3.2.ab_a_a","3.4.ab_c_a","3.4.ab_a_ae"]
dim3_D = ["3.2.ac_c_ac", "3.3.ad_j_ap", "3.2.a_a_ac", "3.7.a_a_abj"]
dim3_E = ["3.3.a_a_aj"]

paper_table_examples = [dim2_A_s, dim2_B_s, dim2_C_s, dim2_A_ns, dim2_B_ns, dim2_C_ns, dim3_A, dim3_B, dim3_C, dim3_D, dim3_E]
paper_table_names = ["dim2_A_s", "dim2_B_s", "dim2_C_s", "dim2_A_ns", "dim2_B_ns", "dim2_C_ns", "dim3_A", "dim3_B", "dim3_C", "dim3_D", "dim3_E"]

# The following function takes as input a list of LMFDB labels and returns a
# string in the format of a magma list, whose entries are the corresponding
# Weil polyinomials
def label_to_poly_Magma_list(name, table_list):
    Pol.<T> = PolynomialRing(QQ)
    return name + " := " + str([weil_poly_from_label(Pol,lab) for lab in table_list]) + ";"

for x in [[paper_table_names[i], paper_table_examples[i]] for i in range(len(paper_table_names))]:
    print(label_to_poly_Magma_list(x[0],x[1]))

