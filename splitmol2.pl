#!/usr/bin/perl
# split a mol2 file into separate files for each residue (one structure/file)

sub next_line()
{
	$nl = <>;
	chomp $nl;
	$nl =~ s/$//; # dos files
	return $nl;
}

# Find atoms header
$line = next_line();
while (($line !~ /^\@<TRIPOS>ATOM/) && (!eof()) ) { $line=next_line(); }
die "<TRIPOS>ATOM header not found.\n" if (eof());

$na = 0;  # nr. of atoms
$nb = 0;  # nr. of bonds
$nr = 0;  # nr. of residues /starts from 1/
@atoms = (); # atom lines
@bonds = (); # residue numbers
@res = ();   # residue numbers of atoms

# read atoms till next header or EOF
$line = next_line();
while (($line !~ /^\@<TRIPOS>/) && (!eof()))
{
	if ($line !~ /^\s*$/) # skip empty lines
	{
		$atoms[$na] = $line;
		@atom = split(' ',$line);
		$res[$na]=$atom[6];
		if ( $res[$na] > $nr ) { $nr = $res[$na]; }
		$na++;
	}
	$line = next_line();
}

# find bonds header
while (($line !~ /^\@<TRIPOS>BOND/) && (!eof()) ) { $line=next_line(); }
# read & store bonds
$line = next_line();
while (($line !~ /^\@<TRIPOS>/) && (!eof()))
{
	if ($line !~ /^\s*$/) # skip empty lines
	{
		# no inter-residue bonds are allowed
		@bond = split(' ',$line);
		die("$ARGV skipped: a bond was found between the residues.\n")
			if ( $res[$bond[1]-1] != $res[$bond[2]-1] );
		$bonds[$nb] = $line;
		$nb++;
	}
	$line = next_line();
}

# find crysinf header - optional
$crys = "";
while (($line !~ /^\@<TRIPOS>CRYSIN/) && (!eof()) ) { $line=next_line(); }
# read & store crysin line
if (!eof()) { $crysin = next_line(); }

# get new atom and bond nubers for the separate files
@num_in_res = (); # nr. of atoms in each residue
@numb_in_res = (); # nr. of bonds in each residue
@newnum = (); # number of the atom within its residue
@newnumb = (); # number of the atom within its residue

for ($i=1; $i<=$nr; $i++)
{
	$n = 0;
	for ($a=0; $a<$na; $a++)
	{
		if ($res[$a] == $i)
		{
			$n++;
			$newnum[$a] = $n;
		}
	}
	$num_in_res[$i] = $n;

	$n = 0;
	for ($b=0; $b<$nb; $b++)
	{
		@bond = split(' ',$bonds[$b]);
		if ($res[$bond[1]-1] == $i)
		{
			$n++;
			$newnumb[$b] = $n;
		}
	}
	$numb_in_res[$i] = $n;
}

# write output files
for ($i=1; $i<=$nr; $i++)
{
	$fn = $ARGV;
	$fn =~ s/\.mol2$//i;
	$mol = "$fn-$i";
	$fn = "$mol.mol2";
	open(OUT,">$fn") || die ("Cannot open file $fn: $!\n");
	print(OUT "@<TRIPOS>MOLECULE\n$mol\n");
	print(OUT "$num_in_res[$i] $numb_in_res[$i] 1\n");
	print(OUT "SMALL\nFORMAL_CHARGES\n****\n");
	print(OUT "created from $ARGV\n\n");

	print(OUT "@<TRIPOS>ATOM\n");
	for ($a=0; $a<$na; $a++)
	{
		if ($res[$a] == $i)
		{
			@atom = split(' ',$atoms[$a]);
			print(OUT "  $newnum[$a] $atom[1] ");
			print(OUT "$atom[2] $atom[3] $atom[4] ");
			print(OUT "$atom[5] 1 RES1 $atom[8]\n");
		}
	}
	print(OUT "@<TRIPOS>BOND\n");
	for ($b=0; $b<$nb; $b++)
	{
		@bond = split(' ',$bonds[$b]);
		if ($res[$bond[1]-1] == $i)
		{
			print(OUT "  $newnumb[$b] ");
			print(OUT "$newnum[$bond[1]-1] $newnum[$bond[2]-1] ");
			print(OUT "$bond[3]\n");
		}
	}

	print(OUT "@<TRIPOS>SUBSTRUCTURE\n    1 RES1 1\n\n");
	if ($crysin ne "")
	{
		print(OUT "@<TRIPOS>CRYSIN\n$crysin\n\n");
	}
	close OUT;
}

