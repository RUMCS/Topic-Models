#------------------------------------------------------------------------------
#This script reads in XML data from StackOverflow and returns it in text format
#Output "title \t body", one line per post
#
#Usage example: perl xml2csv.pl Posts.xml Posts.xml.csv
#------------------------------------------------------------------------------
use warnings;
use strict;

use XML::Simple;    # use XML::Twig for big files
use HTML::Strip;

use utf8;
use Text::Unidecode;    #translates unicode into ascii

my $minVoteCount = 3;   #minimum number of votes required to save an answer

sub transform;

if ( scalar(@ARGV) != 2 ) {
	print "usage: xml2csv.pl in_file_name out_file_name\n";
	exit;
}

#read XML file into memory
my $ref = XMLin( $ARGV[0] );

open( my $fid, ">" . $ARGV[1] )
  or die("Cannot open $ARGV[1] : $!\n");    #export to CSV

my $hs = HTML::Strip->new();
print $fid "create_ts\tid\ttags\ttitle\tbody\tanswers\n";    #header
foreach ( @{ $ref->{row} } ) {
	if ( $_->{PostTypeId} == "1" ) {    #if the post is the original question
		my $creationDate = ( $_->{CreationDate} );
		my $id           = $_->{Id};
		my $tags         = $_->{Tags};
		my $bestAnswerId = $_->{AcceptedAnswerId};

		#remove "<" and replace ">" with " "

		$tags =~ s/></^&^/g;
		$tags =~ s/^<//g;
		$tags =~ s/>$//g;

		my $title   = transform( $_->{Title} );
		my $body    = transform( $hs->parse( $_->{Body} ) ); #strip HTML tags
		my $answers = getAnswersText( $id, $minVoteCount, $bestAnswerId );
		print $fid join( "\t", ( $creationDate, $id, $tags, $title, $body, $answers ) ) . "\n";
	}
}
close $fid;
$hs->eof;

#convert unicode to ascii and make it lower case
sub transform {
	my ($txt) = @_;
	$txt =~ s/\n|\t|\r/ /g;          #replace new line or tab with space
	$txt =~ s/[[:punct:]]/ /g;    #replace punctuation with space
	return lc( unidecode($txt) ); #convert #convert unicode to latin and change words to lowercase
}

sub getAnswersText {
	my ( $questionId, $minVote, $AcceptedAnswerId ) = @_;
	my $answers = "";

	foreach ( @{ $ref->{row} } ) {
		if (
			(
				( $_->{PostTypeId}   eq "2" )            #record type = answer
				and ( $_->{ParentId} eq $questionId )    #answer to question $questionId
				and ( ( $_->{Score} + 0 ) >= $minVote )  #minimal vote
			)
			or (                                         #always get accepted answer
				defined($AcceptedAnswerId) and ( $_->{Id} eq $AcceptedAnswerId )
			)
		  )
		{
			$answers .= transform( $hs->parse( $_->{Body} ) ) . " ";
		}
	}
	return $answers;
}
