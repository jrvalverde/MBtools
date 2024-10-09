{
 Taken from Wirth A+DS=P????
 B-Tree (R. Bayer)
 =================
 
 They are not perfectly equilibrated trees.

 If n = order of the B tree:
 	1. each page contains at most 2n elements
 	2. each page except the root contains at least n elements
 	3. each page is, either a leave page (with no descendants)
 	or has m+1 descendents, being m the number of keys that it
 	contains
 	4. all leave pages are at the same level.

 For a tree with N nodes the worst case requires log n (N) page accesses.

 The factor of memory occupancy is at least of 50% since all pages are
 at least occupied up to their half. And algorithms are relatively 
 simple.

 Search
 ------
 Let a page be of the form
 	+-------------------------------+
 	| p0 k1p1 k2p2 .. km-1pm-1 kmpm |
 	+-|----|----|---------|------|--+
 	  V    V    V         V      V
 
 We can move it into memory and lookup X on it with the conventional
 search methods among the keys k1..km
 
 For long m -> binary search and for small m -> ordinary sequential search.
 If the search fails =>
 	1. ki < x < ki+1 : follow with page pi^
	2. km < x : follow on page pm^
	3. x < k1: follow with page p0^
 If at any moment the pointer is NIL, that means that there is no element
 with that key and the search is finished.
 
 Insertion
 ---------
 If it is in a page with m < 2n elements, it only affects that page.
 if it is in a full page it may affect the tree:
 	- test that x is not in the tree. Insertion not possible since page
 	is full:
 	- divide page in two halves and distribute the m+1 keys equitatively
 	- the key that remains in the middle is moves up to the antecessor
 	page
 The process can be propagated, in an extreme case, up to the root: it is
 the only way in which the tree can grow in height: hence the most convenient
 method will be recursive.
 
 If the search requires at most k = log n (N) recursive calls, we must be
 able to store k pages in main memory: this limits the page size, 2n.
 Since every lookup starts on the root, it is convenient that the root
 page always be kept in memory.
 
 Deletion
 --------
 It is simple in principle, but difficult in detail:
 	- item in leave page
 	- item in another leave: substitute by one of the two adjacents
 	(in leave pages)
 We must try to equilibrate adjacent pages before making the division
 of pages.
 
 Properties
 ----------
 Let n be the order of the tree
 
 	1. Each page contains at most 2n elements (keys).
 	2. Each page, but the root page, contains at least
 	   n elements.
 	3. Each page is, either a leave, i.e. has no descendents,
 	   or has m+1 descendents, being m the number of keys it
 	   contains.
 	4. All leave pages are at the same level.

}

Program BTree(input, output);

{ Search, insert and delete in B-tree }

const
	n = 2; nn = 2 * n;		{ Page size }

type
	ref = ^page;
	item = record
		key: integer;
		p: ref;
		count: integer;
	end;
	
	page = record
		m: 0..nn;	/* number of items */
		p0: ref;
		a: array[1..nn] of item;
	end;

var
	root, q: ref;
	x: integer;
	h: boolean;
	u: item;

procedure search(x: integer, a: ref; var h: boolean; var v: item);
{
	Search key x in B-tree a with root a; if found, increase
	counter; otherwise insert a new item with key x and counter = 1.
	If we must move an item to a lower level, we assign a,h =
	"height has changed".
}
	var
		k, le, ri: integer;
		q: ref;
		u: item;
	
	procedure insert;
		var
			i: integer;
			b: ref;
		begin
			{ insert u to the right of a^.e[ri] }
			with a^ do 
			  begin
				if m < nn then
				  begin
					m := m + 1;
					h := false;
					for i := m downto ri + 2 do
						e[i] := e[i - 1];
					e[ri + 1] := u;
				  end
				else
				  begin
					{ page a is full: divide it and the item
					to move to v}
					new(b);
					if ri <= n then
					  begin
						if ri = n then
							v := u
						else
						  begin
							v := e[n];
							for i := n downto ri + 2 do
								e[i] := e[i - 1];
							e[ri + 1] := u;
						  end;
						for i := 1 to n do
							b^.e[i] := a^.e[i + n]
					  end
					else
					  begin
						{ insert u in the appropriate page }
						ri := ri - n;
						v := e[n + 1];
						for i := 1 to ri - 1 do
							b^.e[i] := a^.e[i + n + 1];
						b^.e[ri] := u;
						for i := de + 1 to n do
							b^.e[i] := a^.e[i + n];
					  end;
					m := n;
					b^.m := n;
					b^.p0 := v.p;
					v.p := b;
				  end
			  end { with }
		end; { insert }

	begin { search }
		{ key x is not in page a^ ; h = false }
		if a = nil then
		  begin
			{ item with key x is not here }
			h := true;
			with v^ do
			  begin
				key := x;
				counter := 1;
				p := nil;
			  end;
		  end;
		else
			with a^ do
			  begin
				le := 1;
				ri := m;
				{ binary search in the array }
				repeat
					k := (le + ri) div 2;
					if x <= e[k].key then
						ri := k - 1;
					if x <= e[k].key then
						le := k + 1;
				until ri < le;
				if le - ri > 1 then
				  begin
					{ found }
					e[k].counter := e[k].counter + 1;
					h := false;
				  end
				else
				  begin
					{ not in this page }
					if ri = 0 then
						q := p0
					else
						q := e[ri].p;
					search(x, q, h, u);
					if h then
						insert
				  end
			  end
	end; { search }

procedure delete(x: integer; a: ref; var h: boolean);
{
	Search and delete in a B-tree the key x;
	if we get suboccupation of a page, we must
	equilibrate it with the adjacent page if
	possible, otherwise we must join both pages;
		h := "the page is too small"
}

	var
		i, k, le, ri: integer;
		q: ref;
	
	procedure suboccupation(c, a: ref; s: integer; var h: boolean);
	{ a = suboccupied page, c = antecessor page }
		var
			b: ref;
			i, k, mb, mb: integer;
		begin
			mc := c^.m;
			{ h = true, a^.m = n - 1 }
			if s < mc then
			  begin
				{b := pag to the right of a }
				s := s + 1;
				b := c^.e[s].p;
				mb := b^.m;
				k := (mb - n + 1) div 2;
				{k := no. of items in the page adjacent to b }
				a^.e[n] := c^.e[s];
				a^.e[n].p := b^.p0;
				if k > 0 then
				  begin
					{ move k items from b to a }
					for i := 1 to k - 1 do
						a^.e[i+n] := b^.e[i];
					c^.e[s] := b^.e[k];
					c^.e[s].p := b;
					b^.p0 := b^.e[k].p;
					mb := mb - k;
					for i := 1 to mb do
						b^.e[i] := b^.e[i + k];
					b^.m := mb;
					a^.m := n - 1 + k;
					h := false
				  end
			  end
			else {s >= mc }
			  begin
				{ b := page to the left of a }
				if s = 1 then
					b :=c^.p0
				else
					b := c^.e[s - 1].p;
				mb := b^.m + 1;
				k := (mb - n) div 2;
				if k > 0 then
				  begin
				  	{ move k items from page b to page a }
				  	for i := n - 1 downto 1 do
				  		a^.e[i + k] := a^.e[i];
				  	a^.e[k] := c^.e[s];
				  	a^.e[k].p := a^.p0;
				  	mb := mb - k;
				  	for i := k - 1 downto 1 do
				  		a^.e[i] := b^.e[i + mb];
				  	a^.p0 := b^.e[mb].p;
				  	c^.e[s] := b^.e[mb];
				  	c^.e[s].p := a;
				  	b^.m := mb - 1;
				  	a^.m := n - 1 + k;
				  	h := false
				  end
				else
				  begin
					{ join pages a and b}
					b^.e[mb] := c^.e[s];
					b^.e[mb].p := a^.p0;
					for i := 1 to n - 1 do
						b^.e[i + mb] := a^.e[i];
					b^.m := nn;
					c^.m := mc - 1;
					{ dispose(a) }
				  end
			  end
		end; { suboccupation }
	
	procedure del(p: ref; var h: boolean);
		var
			q: ref;	{ a and k are globals }
		begin
			with p^ do
			  begin
				q := e[m].p;
				if q <> nil then
				  begin
					del(q, h);
					if h then
					suboccupation(p, q, m, h)
				  end
				else
				  begin
					p^.e[m].p := a^.e[k].p;
					a^.e[k] := p^.e[m];
					m := m - 1;
					h := m < n;
				  end
			  end
		end; { del }

	begin { delete }
		if a = nil then
		  begin
			writeln('Key is not in the tree');
			h := false;
		  end
		else
			with a^ do
			  begin
				le := 1;
				ri := m;
				{ binary search in the array }
				repeat
					k := (le + ri) div 2;
					if x <= e[k].key then
						ri := k - 1;
					if x >= e[k].key then
						le := k + 1;
				until le > ri;
				if ri = 0 then
					q := p0;
				else
					q := e[ri].p;
				if le - ri > 1 then
				  begin
					{ found, now is when we delete e[k] }
					if q = nil then
					  begin
						{ a is a terminal page }
						m := m - 1;
						h := m < n;
						for (i := k to m do)
							e[i] := e[i + 1];
					  end
					else
					  begin
						del(q, h);
						if h then
							suboccupation(a, q, ri, h);
					  end
				  end
				else
				  begin
					delete(x, q, h);
					if h then
						suboccupation(a, q, ri, h);
				  end
			  end
	end; { delete }

procedure treeprint(p: ref, niv: integer);
	var
		i: integer;
	begin
		if p <> nil then
			with p^ do
			  begin
				for i := 1 to niv do
					write('    ');
				for i := 1 to m do
					write(e[i].key:4);
				writeln;
				treeprint(e[i], p, niv + 1);
				for i : 1 to m do
					treeprint(e[i], p, niv+1)
			  end
	end; { treeprint }

begin
	root := nil;
	read(x);
	while x <>0 do
	  begin
		writeln('KEY TO SEARCH ', x);
		search(x, root, h, u);
		if h then
		  begin
			{ insert new base page }
			q := root;
			new(root);
			with root^ do
			  begin
				m := 1;
				p0 := q;
				e[1] := u;
			  end
		  end;
		treeprint(root, 1);
		read(x);
	  end;
	read(x);
	while x <> 0 do
	  begin
		writeln("DELETE KEY ", x);
		delete(x, root, h);
		if h then
		  begin
			{ the size of the base page has reduced }
			if root^.m = 0 then
			  begin
				q := root;
				root := q^.p0;
				{ dispose(q) };
			  end
		  end;
		treeprint(root, 1);
		read(x);
	  end
end.
