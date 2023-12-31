import matplotlib.pyplot as plt

my_id_score_dict = {'131056': (-4, 6), '265933': (4, 5), '587687': (-2, 5), '28413': (6, 6), '289585': (-2, 5), '295165': (-1, 6), '187021': (0, 6), '216031': (-1, 5), '541742': (
    6, 6), '207767': (4, 4), '303454': (5, 6), '425122': (-1, 7), '101503': (5, 6), '220193': (4, 4), '203722': (4, 4), '570826': (-2, 6), '332966': (4, 5), '336665': (4, 4), '250479': (6, 6)}

# my_id_score_dict = {'182109': (-1, 3), '6539': (-3, 6), '245557': (-1, 4), '165418': (5, 5), '408281': (6, 6), '242876': (0, 4), '191799': (5, 5), '310647': (5, 5), '378853': (6, 6), '572612': (-2, 5), '441014': (-1, 6), '259881': (-3, 4), '136380': (-1, 7), '524211': (5, 5), '36320': (-1, 5), '278075': (5, 6), '294254': (-5, 6), '309022': (-1, 6), '78832': (6, 6), '245102': (5, 6), '188127': (5, 5), '262865': (5, 5), '323003': (0, 6), '426536': (5, 7), '274257': (0, 6), '214914': (-3, 3), '301010': (-1, 4), '131056': (-4, 6), '265933': (4, 5), '587687': (-2, 5), '28413': (6, 6), '289585': (-2, 5), '295165': (-1, 6), '187021': (0, 6), '216031': (-1, 5), '541742': (6, 6), '207767': (4, 4), '303454': (5, 6), '425122': (-1, 7), '101503': (5, 6), '220193': (4, 4), '203722': (4, 4), '570826': (-2, 6), '332966': (4, 5), '336665': (4, 4), '250479': (6, 6), '170892': (-1, 6), '304728': (5, 5), '282451': (-1, 4), '104109': (
#     5, 6), '499626': (6, 6), '43917': (-4, 6), '203749': (-1, 4), '16072': (5, 5), '348946': (-1, 4), '510617': (4, 4), '563236': (6, 6), '108318': (5, 5), '514212': (-1, 5), '140216': (-1, 6), '156444': (6, 6), '183577': (6, 6), '248112': (4, 4), '549034': (-1, 3), '9459': (5, 6), '242651': (-2, 5), '439022': (6, 6), '579169': (-1, 4), '94867': (-1, 7), '570233': (-1, 7), '353130': (4, 5), '509956': (-1, 2), '191389': (4, 6), '520891': (4, 4), '275907': (5, 6), '25577': (-1, 6), '135405': (4, 5), '252737': (-2, 6), '279026': (-1, 6), '356995': (6, 6), '445533': (-6, 6), '544153': (6, 7), '328359': (5, 6), '529772': (-1, 6), '194230': (-1, 5), '240360': (-1, 3), '163571': (-1, 6), '296330': (-1, 5), '188792': (5, 6), '145434': (5, 5), '521447': (5, 5), '272289': (5, 5), '270303': (5, 5), '108379': (-2, 6), '216520': (4, 6), '238495': (-1, 5), '275661': (-1, 5), '109889': (-1, 2), '98692': (-1, 6), '316765': (-1, 6)}

keys = list(my_id_score_dict.keys())
values = list(my_id_score_dict.values())

x_pos = list(range(len(keys)))
# x_pos = keys

# plt.bar([i+0.2 for i in x_pos], [v[0]
#         for v in values], width=0.4, color='powderblue', align='center')

plt.bar([i-0.2 for i in x_pos], [v[1] if v[1] >= 0 else 0 for v in values],
        width=0.3, color='dimgray', align='center', label='True known taxa')

plt.bar([i+0.2 for i in x_pos], [abs(v[0]) if v[0] <
        0 else 0 for v in values], width=0.3, color='maroon', align='center', label='Guessed taxa, with mismatch')

plt.bar([i+0.2 for i in x_pos], [v[0] if v[0] >=
        0 else 0 for v in values], width=0.3, color='darkolivegreen', align='center', label='Guessed taxa, no mismatch')

# plt.bar([i-0.2 for i in x_pos], [v[1]
#         for v in values], width=0.4, color='peachpuff', align='center')

plt.xlabel("ID")
plt.ylabel("Rank (taxon)")
plt.title("ID Scores")


print(keys)
print(x_pos)
plt.xticks(x_pos, keys, rotation=90)
# plt.yticks(range(8))
plt.yticks(range(8), ["", "kingdom", "phylum", "class",
           "order", "family", "genus", "species"])
plt.legend()

plt.show()
