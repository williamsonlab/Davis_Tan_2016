#!/usr/bin/perl

use Tree::Simple;

#
# modify_vX.pl
#
# ./modify_vX.pl --rna=Z --carbamido=X --modfile_specific=STRING --modfile_wildcard=STRING 
#
# rna 0 - assume protein (DEFAULT)
# rna 1 - assume RNA, so ignore the carbamido 
#
# carbamido 0 - ignore carbamido modified cysteines
# carbamido 1 - include only carbamido modified cysteines (DEFAULT)
# carbamido 2 - include carbamido unmodified and modified cysteines
#
# MODFILE_SPECIFIC
# Protein AA Mod Include
# S5      1  a   1
# S12     88 b   1
# S12     88 d   0
# S16     23 m   2
#
# Flag
# 0 - don't include the modification
# 1 - always include the modification
# 2 - include both modified and unmodified versions
#
# MODFILE_WILDCARD
# AA	Mod		Qualifier	Include
# C		c		*			1			-> equivalent to standard carbadmido cysteines
# K		m		*			2			-> methylate(m) every(*) lysine(K) both versions(2)
# H		m		N			1			-> methylate N-terminal hystidines always
# W		m		C			1			-> methylate C-terminal tryptophan always
# *		m		N			1			-> methylte any N-terminal AA always
#
# OUTPUT
# protein seq startres endres missed mod\n";
#

use Getopt::Long;
require "$ENV{'MASSACRE_PATH'}/msc/mikesubs.pl";

$rna = 0;
$carbamido = 1;
$modfile_specific = "";
$modfile_wildcard = "";


&GetOptions(
	"rna=i" => \$rna,
	"carbamido=i" => \$carbamido,
	"modfile_specific=s" => \$modfile_specific,
	"modfile_wildcard=s" => \$modfile_wildcard,
);


print "#RNA $rna\n";
print "#CARBAMIDO $carbamido\n";
print "#MODFILE_SPECIFIC $modfile_specific\n";
print "#MODFILE_WILDCARD $modfile_wildcard\n";

#
# $mod_specific{$protein}[$aa_number][0] -> array of modifications
#	-each array element is the actual modification (m, c, whatever)
# $mod_specific{$protein}[$aa_number][1] -> array of flags
#	-each array element is the flag (0:never,1:always,2:sometimes)
#
if($modfile_specific ne "")
{
	@modfile_specific_array = readfile($modfile_specific);
	@array = split ' ', $modfile_specific_array[0];
	
	if(lc($array[0]) eq "protein")
	{
		$header = shift(@modfile_specific_array);
	}
	
	for($ctr = 0; $ctr <= $#modfile_specific_array; $ctr++)
	{
		@array = split ' ', $modfile_specific_array[$ctr];

		$protein = shift(@array);
		$aa_number = shift(@array);
		$mod = shift(@array);
		$flag = shift(@array);
		
		# Hash -> Array -> Array -> Array
		push(@{$mod_specific{$protein}[$aa_number][0]}, $mod);
		push(@{$mod_specific{$protein}[$aa_number][1]}, $flag);
	}
}

#
# $mod_wildcard{$aa}[0] -> array of modifications
#	-each array element is the actual modification (m, c, whatever)
# $mod_wildcard{$aa}[1] -> array of flags
#	-each array element is the flag (0:never,1:always,2:sometimes)
# $mod_wildcard{$aa}[2] -> array of qualifiers
#	-each array element is a qualifier applied to the aa (*:all,N:nterm,C:cterm)
#
if($modfile_wildcard ne "")
{
	@modfile_wildcard_array = readfile($modfile_wildcard);
	@array = split ' ', $modfile_wildcard_array[0];
	
	if(lc($array[0]) eq "aa")
	{
		$header = shift(@modfile_wildcard_array);
	}
	
	for($ctr = 0; $ctr <= $#modfile_wildcard_array; $ctr++)
	{
		@array = split ' ', $modfile_wildcard_array[$ctr];

		$aa = shift(@array);
		$mod = shift(@array);
		$qualifier = shift(@array);
		$flag = shift(@array);
		
		push(@{$mod_wildcard{$aa}[0]}, $mod);
		push(@{$mod_wildcard{$aa}[1]}, $flag);
		push(@{$mod_wildcard{$aa}[2]}, $qualifier);
	}
}

# This just folds the standard carbamido modification of cysteine into the wildcard system
if($carbamido == 1 && $rna == 0)
{
	push(@{$mod_wildcard{"C"}[0]}, "c");
	push(@{$mod_wildcard{"C"}[1]}, 1);
	push(@{$mod_wildcard{"C"}[2]}, "*");
}
elsif($carbamido == 2 && $rna == 0)
{
	push(@{$mod_wildcard{"C"}[0]}, "c");
	push(@{$mod_wildcard{"C"}[1]}, 2);
	push(@{$mod_wildcard{"C"}[2]}, "*");
}


#
# Read in the output from enzyme_vX.pl
# protein seq startres endres missed
#
$header_flag = 0;
while(defined($input=<STDIN>))
{
	chomp($input);
	
	if(substr($input, 0, 1) eq "#")
	{
		print "$input\n";
		next;
	}

	if($header_flag == 0)
	{
		print "protein seq startres endres missed mod\n";
		$header_flag = 1;
	}
	
	@array = split ' ', $input;
	if(lc($array[0]) eq "protein")
	{
		$header = $input;
		next;
	}
	
	$protein = shift(@array);
	$seq = shift(@array);
	$startres = shift(@array);
	$endres = shift(@array);
	$missed = shift(@array);
	
	@seq_array = split '', $seq;
	$pep_length = $#seq_array + 1;
	
	@mod_array = ();
	@flag_array = ();
	
	@mod_array2 = ();
	$mod_ctr = 0;
	
	
	#
	# Start with the last residue as the root and build the tree towards the first residue
	# ACGK -> K(root) -> G -> C -> A
	# The printing at the end is done starting at the leaf (A), so this means that the final output is in the correct order
	#
	#
	#
	my $tree = Tree::Simple->new("ROOT", Tree::Simple->ROOT);
	
	@old_nodes = ($tree);
	
	# First Pass ($ctr == $endres)
	#	@old_nodes = ($tree)
	#	Add new nodes to the root node ($tree)
	#	These nodes become @new_nodes
	#	At the end @old_nodes is set to @new_nodes.
	
	for($ctr = $endres; $ctr >= $startres; $ctr--)
	{
		@peptide_mods = ();
		
		$residue = pop(@seq_array);
		
		# Specific Modifications
		for($ctr2 = 0; $ctr2 <= $#{$mod_specific{$protein}[$ctr][0]}; $ctr2++)
		{
			if($mod_specific{$protein}[$ctr][1][$ctr2] == 0)
			{
				# Skip this one
				next;
			}
			elsif($mod_specific{$protein}[$ctr][1][$ctr2] == 1)
			{
				# Actual modification
				push(@{$peptide_mods[$#peptide_mods+1]}, $mod_specific{$protein}[$ctr][0][$ctr2]);
				
				@new_nodes = ();
				for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
				{
					$node = Tree::Simple->new($mod_specific{$protein}[$ctr][0][$ctr2]);
					$old_nodes[$ctr3]->addChild($node);
					push(@new_nodes, $node);
				}
				@old_nodes = @new_nodes;
				
			}
			elsif($mod_specific{$protein}[$ctr][1][$ctr2] == 2)
			{
				# Null modification
				push(@{$peptide_mods[$#peptide_mods+1]}, "#");
				# Actual modification
				push(@{$peptide_mods[$#peptide_mods+1]}, $mod_specific{$protein}[$ctr][0][$ctr2]);
				
				@new_nodes = ();
				for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
				{
					$node = Tree::Simple->new("#");
					$old_nodes[$ctr3]->addChild($node);
					push(@new_nodes, $node);
					$node = Tree::Simple->new($mod_specific{$protein}[$ctr][0][$ctr2]);
					$old_nodes[$ctr3]->addChild($node);
					push(@new_nodes, $node);
				}
				@old_nodes = @new_nodes;
			}
		}
		
		# Wildcard modifications (Specific Residue)
		for($ctr2 = 0; $ctr2 <= $#{$mod_wildcard{$residue}[0]}; $ctr2++)
		{
			if($mod_wildcard{$residue}[1][$ctr2] == 0)
			{
				# Skip this one
				next;
			}
			elsif($mod_wildcard{$residue}[1][$ctr2] == 1)
			{
				if($mod_wildcard{$residue}[2][$ctr2] eq "*")
				{
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
				elsif($mod_wildcard{$residue}[2][$ctr2] eq "N" && $ctr == $startres)
				{
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
				elsif($mod_wildcard{$residue}[2][$ctr2] eq "C" && $ctr == $endres)
				{
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
			}
			elsif($mod_wildcard{$residue}[1][$ctr2] == 2)
			{
				if($mod_wildcard{$residue}[2][$ctr2] eq "*")
				{
					# Null modification
					push(@{$peptide_mods[$#peptide_mods+1]}, "#");	
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);	

					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new("#");
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
				elsif($mod_wildcard{$residue}[2][$ctr2] eq "N" && $ctr == $startres)
				{
					# Null modification
					push(@{$peptide_mods[$#peptide_mods+1]}, "#");	
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);	

					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new("#");
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
				elsif($mod_wildcard{$residue}[2][$ctr2] eq "C" && $ctr == $endres)
				{
					# Null modification
					push(@{$peptide_mods[$#peptide_mods+1]}, "#");	
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);	
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new("#");
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
			}
		}
		
		# Wildcard modifications (Wildcard Residue)
		$residue = "*";
		for($ctr2 = 0; $ctr2 <= $#{$mod_wildcard{$residue}[0]}; $ctr2++)
		{
			if($mod_wildcard{$residue}[1][$ctr2] == 0)
			{
				# Skip this one
				next;
			}
			elsif($mod_wildcard{$residue}[1][$ctr2] == 1)
			{
				if($mod_wildcard{$residue}[2][$ctr2] eq "*")
				{
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
				elsif($mod_wildcard{$residue}[2][$ctr2] eq "N" && $ctr == $startres)
				{
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
				elsif($mod_wildcard{$residue}[2][$ctr2] eq "C" && $ctr == $endres)
				{
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
			}
			elsif($mod_wildcard{$residue}[1][$ctr2] == 2)
			{
				if($mod_wildcard{$residue}[2][$ctr2] eq "*")
				{
					# Null modification
					push(@{$peptide_mods[$#peptide_mods+1]}, "#");	
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);	
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new("#");
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
				elsif($mod_wildcard{$residue}[2][$ctr2] eq "N" && $ctr == $startres)
				{
					# Null modification
					push(@{$peptide_mods[$#peptide_mods+1]}, "#");	
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);	
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new("#");
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
				elsif($mod_wildcard{$residue}[2][$ctr2] eq "C" && $ctr == $endres)
				{
					# Null modification
					push(@{$peptide_mods[$#peptide_mods+1]}, "#");	
					# Actual modification
					push(@{$peptide_mods[$#peptide_mods+1]}, $mod_wildcard{$residue}[0][$ctr2]);	
					
					@new_nodes = ();
					for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
					{
						$node = Tree::Simple->new("#");
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
						$node = Tree::Simple->new($mod_wildcard{$residue}[0][$ctr2]);
						$old_nodes[$ctr3]->addChild($node);
						push(@new_nodes, $node);
					}
					@old_nodes = @new_nodes;
				}
			}
		}
		
		#
		# Inter-residue delimiter
		#
		# We also trigger this condition if the peptide is of length 1 ("R" or something like that) and there are no modifications
		# No modifications and 1 length means a null tree and no output
		# Really, we don't care about such peptides but it's necessary for completeness
		#
		if($ctr > $startres || $#peptide_mods == -1)
		{				
			push(@{$peptide_mods[$#peptide_mods+1]}, "-");
			
			@new_nodes = ();
			for($ctr3 = 0; $ctr3 <= $#old_nodes; $ctr3++)
			{
				$node = Tree::Simple->new("-");
				$old_nodes[$ctr3]->addChild($node);
				push(@new_nodes, $node);
			}
			@old_nodes = @new_nodes;
		}		
	}
	
	#
	# Output
	#
	$tree->traverse(sub {
		my ($_tree) = @_;
				
		if($_tree->isLeaf())
		{			
			$leaf_output = $_tree->getNodeValue();
			#print "", $_tree->getNodeValue();
			print_parent($_tree);
			
			#print "\n";
			$leaf_output =~ s/\#//g; # Strip out null modifications
			
			if($leaf_output eq "")
			{
				$pep_length = $startres - $endres + 1;
				
				if($pep_length == 1)
				{
					$leaf_output = "-";
				}
				else
				{
					die "Error, incorrect modification string |$leaf_output| for $seq\n";
				}
			}
			
			print "$input $leaf_output\n";
		}
	});
}

sub print_parent
{
	my ($_tree) = @_;
	$parent = $_tree->getParent();
	if($parent->isRoot)
	{
		return;
	}

	$leaf_output = join '', $leaf_output, $parent->getNodeValue();
	#print "", $parent->getNodeValue();
	
	print_parent($parent);
}
