###############################################
#########Python code for range paper###########


cd ~/Desktop
conda activate myenv
python3
import numpy as np
import msprime, pyslim
import tskit
import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd



#Load tree sequences, recapitate, add mutations

ts = pyslim.load("no_contraction.trees")
recap_ts = ts.recapitate(recombination_rate=1e-9, Ne=4000)
ts = pyslim.SlimTreeSequence(msprime.mutate(recap_ts, rate=1e-8))
ts.dump("no_contraction_mut.trees")
ts = pyslim.load("no_contraction_mut.trees")
	

#################################
####Methods for getting pi######

#Calculate heterozygosity across space and time

#times = [800, 650, 550, 450, 350, 0]
#Below is a lot of repetitive code meant to make files for each time-point, including ind. locs, het, and relative het; these are joined later in R

#time 800
alive = ts.individuals_alive_at(800)
targets = np.random.choice(alive, size=100, replace=False)
target_node_list = []
new_nodes = []
k = 0
for ind in targets:
    target_node_list.extend(ts.individual(ind).nodes)
    new_nodes.append([k, k+1])
    k += 2

sts = ts.simplify(target_node_list)
new_targets = [sts.node(u[0]).individual for u in new_nodes]
het = sts.diversity([sts.individual(u).nodes for u in new_targets])
het_df = pd.DataFrame(het, columns = ["het"])
relative = het / np.mean(het)
relative_df = pd.DataFrame(relative, columns = ["relative"])

#write locs and het to a file

indivlist = []
indivnames = []
with open("frag_800.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y"]) + "\n")
  for i in new_targets:
        indivlist.append(i)
        ind = sts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label,
                str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")

df = pd.read_csv('frag_800.txt', delimiter = "\t")
df["het"] = het_df
df["time"] = 800
df["relative"] = relative_df
df.to_csv("frag_800.txt", index=False)

#time 650
alive = ts.individuals_alive_at(650)
targets = np.random.choice(alive, size=100, replace=False)
target_node_list = []
new_nodes = []
k = 0
for ind in targets:
    target_node_list.extend(ts.individual(ind).nodes)
    new_nodes.append([k, k+1])
    k += 2

sts = ts.simplify(target_node_list)
new_targets = [sts.node(u[0]).individual for u in new_nodes]
het1 = sts.diversity([sts.individual(u).nodes for u in new_targets])
het_df = pd.DataFrame(het1, columns = ["het"])
relative = het1 / np.mean(het)
relative_df = pd.DataFrame(relative, columns = ["relative"])

#write locs and het to a file

indivlist = []
indivnames = []
with open("frag_650.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y"]) + "\n")
  for i in new_targets:
        indivlist.append(i)
        ind = sts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label,
                str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")

df = pd.read_csv('frag_650.txt', delimiter = "\t")
df["het"] = het_df
df["time"] = 650
df["relative"] = relative_df
df.to_csv("frag_650.txt", index=False)

#time 550
alive = ts.individuals_alive_at(550)
targets = np.random.choice(alive, size=100, replace=False)
target_node_list = []
new_nodes = []
k = 0
for ind in targets:
    target_node_list.extend(ts.individual(ind).nodes)
    new_nodes.append([k, k+1])
    k += 2

sts = ts.simplify(target_node_list)
new_targets = [sts.node(u[0]).individual for u in new_nodes]
het1 = sts.diversity([sts.individual(u).nodes for u in new_targets])
het_df = pd.DataFrame(het1, columns = ["het"])
relative = het1 / np.mean(het)
relative_df = pd.DataFrame(relative, columns = ["relative"])

#write locs and het to a file

indivlist = []
indivnames = []
with open("frag_550.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y"]) + "\n")
  for i in new_targets:
        indivlist.append(i)
        ind = sts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label,
                str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")

df = pd.read_csv('frag_550.txt', delimiter = "\t")
df["het"] = het_df
df["time"] = 550
df["relative"] = relative_df
df.to_csv("frag_550.txt", index=False)

#time 450
alive = ts.individuals_alive_at(450)
targets = np.random.choice(alive, size=100, replace=False)
target_node_list = []
new_nodes = []
k = 0
for ind in targets:
    target_node_list.extend(ts.individual(ind).nodes)
    new_nodes.append([k, k+1])
    k += 2

sts = ts.simplify(target_node_list)
new_targets = [sts.node(u[0]).individual for u in new_nodes]
het1 = sts.diversity([sts.individual(u).nodes for u in new_targets])
het_df = pd.DataFrame(het1, columns = ["het"])
relative = het1 / np.mean(het)
relative_df = pd.DataFrame(relative, columns = ["relative"])

#write locs and het to a file

indivlist = []
indivnames = []
with open("frag_450.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y"]) + "\n")
  for i in new_targets:
        indivlist.append(i)
        ind = sts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label,
                str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")

df = pd.read_csv('frag_450.txt', delimiter = "\t")
df["het"] = het_df
df["time"] = 450
df["relative"] = relative_df
df.to_csv("frag_450.txt", index=False)

#time 350
alive = ts.individuals_alive_at(350)
targets = np.random.choice(alive, size=100, replace=False)
target_node_list = []
new_nodes = []
k = 0
for ind in targets:
    target_node_list.extend(ts.individual(ind).nodes)
    new_nodes.append([k, k+1])
    k += 2

sts = ts.simplify(target_node_list)
new_targets = [sts.node(u[0]).individual for u in new_nodes]
het1 = sts.diversity([sts.individual(u).nodes for u in new_targets])
het_df = pd.DataFrame(het1, columns = ["het"])
relative = het1 / np.mean(het)
relative_df = pd.DataFrame(relative, columns = ["relative"])

#write locs and het to a file

indivlist = []
indivnames = []
with open("frag_350.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y"]) + "\n")
  for i in new_targets:
        indivlist.append(i)
        ind = sts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label,
                str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")

df = pd.read_csv('frag_350.txt', delimiter = "\t")
df["het"] = het_df
df["time"] = 350
df["relative"] = relative_df
df.to_csv("frag_350.txt", index=False)

#time 0
alive = ts.individuals_alive_at(0)
targets = np.random.choice(alive, size=100, replace=False)
target_node_list = []
new_nodes = []
k = 0
for ind in targets:
    target_node_list.extend(ts.individual(ind).nodes)
    new_nodes.append([k, k+1])
    k += 2

sts = ts.simplify(target_node_list)
new_targets = [sts.node(u[0]).individual for u in new_nodes]
het1 = sts.diversity([sts.individual(u).nodes for u in new_targets])
het_df = pd.DataFrame(het1, columns = ["het"])
relative = het1 / np.mean(het)
relative_df = pd.DataFrame(relative, columns = ["relative"])

#write locs and het to a file

indivlist = []
indivnames = []
with open("frag_0.txt", "w") as indfile:
  indfile.writelines("\t".join(["vcf_label"]
                               + ["x", "y"]) + "\n")
  for i in new_targets:
        indivlist.append(i)
        ind = sts.individual(i)
        vcf_label = f"tsk_{ind.id}"
        indivnames.append(vcf_label)
        data = [vcf_label,
                str(ind.location[0]), str(ind.location[1])]
        indfile.writelines("\t".join(data) + "\n")

df = pd.read_csv('frag_0.txt', delimiter = "\t")
df["het"] = het_df
df["time"] = 0
df["relative"] = relative_df
df.to_csv("frag_0.txt", index=False)



###############################
#####Plotting IBD########


np.random.seed(111)

alive = ts.individuals_alive_at(0)
locs = ts.individual_locations[alive, :]
old_ones = ts.individuals_alive_at(800)
old_locs = ts.individual_locations[old_ones, :]

#Check where folks are

xmax = max(locs[:,0])
ymax = max(locs[:,1])
xmin = min(locs[:, 0])
ymin = min(locs[:,1])

L = ymax - ymin
w = xmax - xmin


#shrinkage groups
W = 20
w = 3
groups = {
   'bottomleft' : alive[np.logical_and(locs[:, 0] < 9, locs[:, 1] < 8)],
   'topleft' : alive[np.logical_and(locs[:, 0] < 8, locs[:, 1] > 12)],
   'bottomright' : alive[np.logical_and(locs[:, 0] > 10, locs[:, 1] < 8)],
   'topright' : alive[np.logical_and(locs[:, 0] > 10, locs[:, 1] > 12)],
   'center' : alive[np.logical_and(np.abs(locs[:, 0] - W/2) < w/2,
                                   np.abs(locs[:, 1] - W/2) < w/2)],
   'ancient' : np.random.choice(old_ones, size=50)}

#amputation groups
groups = {
   'top' : alive[np.logical_and(locs[:, 0] < 3, locs[:,1] < ymax)],
   'upper_middle' : alive[np.logical_and(locs[:,0] < 7, locs[:,0] > 4)],
   'middle' : alive[np.logical_and(locs[:,0] < 11, locs[:,0] > 8)],
   'lower_middle' : alive[np.logical_and(locs[:,0] > 13, locs[:,0] < 16)],
   'lower' : alive[np.logical_and(locs[:,0] < xmax, locs[:,0] > 17)],
   'ancient' : np.random.choice(old_ones, size=50)}

#frag groups
groups = {
   'bottomleft' : alive[np.logical_and(locs[:, 0] < w, locs[:, 1] < w)],
   'topleft' : alive[np.logical_and(locs[:, 0] < 7.5, locs[:, 1] > W - w)],
   'bottomright' : alive[np.logical_and(locs[:, 0] > 15, locs[:, 1] < 5)],
   'topright' : alive[np.logical_and(locs[:, 0] > 18, locs[:, 1] > 16)],
   'ancient' : np.random.choice(old_ones, size=50)}

for k in groups:
	print(f"We have {len(groups[k])} individuals in the {k} group")


#shrinkage group order
group_order = ['bottomleft', 'topleft', 'bottomright', 'topright', 'center', 'ancient']
#amputation group order
group_order = ['top', 'upper_middle', 'middle', 'lower_middle', 'lower', 'ancient']
#frag group order
group_order = ['topleft', 'topright', 'bottomleft', 'bottomright', 'ancient']

ind_colors = np.repeat(0, ts.num_individuals)
for j, k in enumerate(group_order):
   ind_colors[groups[k]] = 1 + j

old_locs = ts.individual_locations[old_ones, :]

#Plot where sampled individuals are

fig = plt.figure(figsize=(12, 6), dpi=300)
ax = fig.add_subplot(121)
ax.set_title("Timestep 0")
ax.scatter(locs[:,0], locs[:,1], s=20, c=ind_colors[alive])
ax.set_xlim(0,20)
ax.set_ylim(0,20)
ax = fig.add_subplot(122)
ax.set_title("Timestep 800")
ax.scatter(old_locs[:, 0], old_locs[:, 1], s=20, c=ind_colors[old_ones])
ax.set_xlim(0,20)
ax.set_ylim(0,20)
fig.savefig("sampled_locs.png")

sampled_nodes = [[] for _ in groups]
for j, k in enumerate(group_order):
   for ind in groups[k]:
      sampled_nodes[j].extend(ts.individual(ind).nodes)

pairs = [(i, j) for i in range(5) for j in range(5)]
group_div = ts.divergence(sampled_nodes, indexes=pairs).reshape((5, 5))

print("\t" + "\t".join(group_order))
for i, group in enumerate(group_order):
   print(f"{group_order[i]}:\t" + "\t".join(map(str, np.round(group_div[i], 7))))

Fst = ts.Fst(sampled_nodes, indexes = pairs).reshape((5, 5))

print("\t" + "\t".join(group_order))
for i, group in enumerate(group_order):
   print(f"{group_order[i]}:\t" + "\t".join(map(str, np.round(Fst[i], 4))))


##############################
####Plot ancestry spread######


sys.path.append("/Users/zacharyhancock/Desktop/spgr/plots")
import spatial_slim as sps
import plot_ancestry_spread
import scipy
import matplotlib._pylab_helpers

def plot_ancestry(ts, targets, times, node_ancestry):
    locs = ts.individual_locations
    def size_fun(inds, k, scale=2000):
        return np.fromiter(map(lambda x: scale * sum(node_ancestry[k][ts.individual(x).nodes]), inds),
                           'float')
    xmax = max(locs[:,0])
    ymax = max(locs[:,1])
    figs = []
    for time in times:
        fig = plt.figure(figsize=(6, 6 * 20 / 20))
        ax = fig.add_subplot(111)
        plt.axis('equal')
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 20)
        sps.plot_density(ts, time, ax, scatter=False)
        colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a","#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#E76060","#1F647F"]
        inds = ts.individuals_by_time(time)
        for k, child in enumerate(targets):
            ax.scatter(locs[inds, 0],
                       locs[inds, 1], 
                       sizes=size_fun(inds, k),
                       facecolor=colors[k],
                       edgecolor=None,
                       alpha=0.75)
        figs.append(fig)
    return figs

ts = sps.SpatialSlimTreeSequence(pyslim.load("no_contraction_mut.trees"), dim=2)

times = [350, 800]
locs = ts.individual_locations
alive = ts.individuals_alive_at(350)
targets = list(np.random.choice(np.where(alive)[0], 4))
children_nodes = [ts.individual(ind).nodes for ind in targets]
node_ancestry = ts.proportion_ancestry_nodes(children_nodes)

figs = plot_ancestry(ts, targets, times, node_ancestry)

figures=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]

for i, figure in enumerate(figures):
    figure.savefig('no_contraction_figure%d.png' % i)

exit()


#TMRCA and extent of LD

ts = pyslim.load("amputation_mut.trees")

alive = ts.individuals_alive_at(450)
targets = np.random.choice(alive, size=100, replace=False)
target_node_list = []
new_nodes = []
k = 0
for ind in targets:
    target_node_list.extend(ts.individual(ind).nodes)
    new_nodes.append([k, k+1])
    k += 2

sts = ts.simplify(target_node_list)
tmrca = np.zeros(sts.num_trees)
breakpoints = np.zeros(sts.num_trees)
for tree in sts.trees():
    tmrca[tree.index] = tree.time(tree.root)
    breakpoints[tree.index] = tree.interval[0]

tmrca_df = pd.DataFrame(data = tmrca, columns = ['TMRCA'])
breakpoint_df = pd.DataFrame(data = breakpoints, columns = ['breakpoint'])
combined_df = pd.concat([tmrca_df, breakpoint_df], axis=1)
combined_df.to_csv('no_contraction_tmrca_during.csv', sep=' ', index=True)

