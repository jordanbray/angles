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
	sed 's/2^(3\/2)/pt.clone()/g' | \
	sed 's/2^(7\/2)/ps.clone()/g' | \
	sed 's/2^(9\/2)/pq.clone()/g' | \
	sed 's/2^(17\/2)/pr.clone()/g' | \
	sed 's/2^(19\/2)/pu.clone()/g' | \
	sed 's/2^(25\/2)/pv.clone()/g' | \
	sed 's/2^(41\/2)/pw.clone()/g' | \
	sed 's/2^(37\/2)/px.clone()/g' | \
	sed 's/2^(27\/2)/py.clone()/g' | \
	sed 's/2^(21\/2)/pz.clone()/g' | \
	sed -r 's/[0-9]+/nu("\0")/g' | \
	sed 's/^/    let /' | \
	sed 's/$/;/'

