cat test.maxima | maxima | \
	sed 's/(%o[0-9]*) *//g' | \
	sed '/(%i[0-9]*)/d' | \
	sed 's/sqrt(2)/sq.clone()/g' | \
	sed 's/sqrt(sq.clone()+2)/p.clone()/g' | \
	sed 's/sqrt(2+sq.clone())/p.clone()/g' | \
	sed 's/sqrt(2-sq.clone())/n.clone()/g' | \
	sed 's/sqrt(2-n.clone())/nn.clone()/g' | \
	sed 's/sqrt(p.clone()+2)/pp.clone()/g' | \
	sed 's/sqrt(2+p.clone())/pp.clone()/g' | \
	sed 's/sqrt(2-p.clone())/np.clone()/g' | \
	sed 's/sqrt(n.clone()+2)/pn.clone()/g' | \
	sed 's/sqrt(2+n.clone())/pn.clone()/g' | \
	sed 's/c_a/c_a.clone()/' | \
	sed 's/c_b/c_b.clone()/' | \
	sed 's/c_c/c_c.clone()/' | \
	sed 's/c_d/c_d.clone()/' | \
	sed 's/c_e/c_e.clone()/' | \
	sed 's/c_f/c_f.clone()/' | \
	sed 's/c_g/c_g.clone()/' | \
	sed 's/c_h/c_h.clone()/' | \
	sed -re 's/2\^[(]([0-9]+)\/2[)]/twos(\1)/g' | \
	sed -r 's/[0-9]+/nu("\0")/g' | \
	sed -r 's/twos[(]nu[(]"([0-9]+)"[)][)]/twos(\1)/g' | \
	sed 's/^/    let /' | \
	sed 's/$/;/'

