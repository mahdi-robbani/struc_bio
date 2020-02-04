fetch 2ptc
hide everything, resn hoh
hide cartoon
show line
select SC, sidechain
hide line, SC

color skyblue, chain E
color darksalmon, chain I

select active_site, (sidechain or name ca) and (chain E and resi 57+102+195)
show stick, active_site
color cyan, active_site

select joined_inhib, chain I and resi 14+15+16+36+37
show stick, joined_inhib
color ruby, joined_inhib

select inhibitor, chain I and (active_site around 4)
color red, inhibitor

label (active_site and name ca), "%s%s" % (resn,resi)
set label_color, black, (active_site and name ca)
label (joined_inhib and name ca), "%s%s" % (resn,resi)
set label_color, black, (joined_inhib and name ca)

zoom active_site, 4

ray 650,650
png handin1.png