#!/usr/bin/ruby

# require 'minitest/unit'
require 'minitest/autorun'
require './GFF3Parser.rb'

class TC_GFF3Parser < MiniTest::Unit::TestCase

  class Test_GFF3Parser
    include GFF3Parser
  end

  def setup
    @path2speAscafa = "sample_data/spur_scaffold.fas"
    @path2speAgff = "./sample_data/spur_sample.gff"
    @testgenehash = { "gene" => [] }
    @testtranscripthash = { "CDS"                      => [],
                            "exon"                     => [],
                            "intron"                   => [],
                            "start_codon"              => [],
                            "stop_codon"               => [],
                            "mRNA"                     => [],
                            "transcription_end_site"   => [],
                            "transcription_start_site" => [],
                            "3-utr"                    => [],
                            "5-utr"                    => [] }
    @g3p = Test_GFF3Parser::GFF3Parsed.new( @path2speAgff, @path2speAscafa )
  end

  def test_sampletest
    assert_equal( "good result", @g3p.sampletest )
  end

  def test_print_usage
    assert_equal( "\n ./GFF3Parser.rb --gff      <gff3> \n \
                        --scaffold <project_id> \n\n ", @g3p.print_usage )
  end

  def test_create_gene_hash
    assert_equal( @testgenehash, @g3p.create_gene_hash )
  end

  def test_create_transcript_hash
    assert_equal( @testtranscripthash, @g3p.create_transcript_hash )
  end

  def test_parse_gff3
    ftype = "gene"
    line  = "Scaffold881     GLEAN3-UTR      gene    59053   94582   .       -       .       ID=SPU_025302gn;Name=SPU_025302;"
    h = {}
    assert_equal( ["SPU_025302gn"], @g3p.parse_gff3( line, ftype, h ).keys )
    assert_equal( "Scaffold881",    @g3p.parse_gff3( line, ftype, h )["SPU_025302gn"]["gene"][0] )
    assert_equal( "59053", @g3p.parse_gff3( line, ftype, h )["SPU_025302gn"]["gene"][3] )
    assert_equal( "94582", @g3p.parse_gff3( line, ftype, h )["SPU_025302gn"]["gene"][4] )
    assert_equal( "-",     @g3p.parse_gff3( line, ftype, h )["SPU_025302gn"]["gene"][6] )
  end

  def test_create_scaffold_hash
    assert_equal( 397093, @g3p.scaffoldh["Scaffold881"].size )
  end

  def test_create_gff_hash
    assert_equal( ["gene", "SPU_025302-tr"], @g3p.gffh["SPU_025302gn"].keys )
  end

  def test_parse_exon_sites
    @g3p.agene( "SPU_025302gn" ).atranscript( "SPU_025302-tr" ).sort_exon_sites
    exon_values = @g3p.agene( "SPU_025302gn" ).atranscript( "SPU_025302-tr" ).parsed_exon_sites
    assert_equal( [[61717, 62175], [63168, 63380], [64349, 64480], [88680, 88763], [89086, 89453], [94546, 94582]], exon_values  )
  end

  def test_parse_cds_sites
    @g3p.agene( "SPU_025302gn" ).atranscript( "SPU_025302-tr" ).sort_cds_sites
    cds_values  = @g3p.agene( "SPU_025302gn" ).atranscript( "SPU_025302-tr" ).parsed_cds_sites
    assert_equal( [61717, 62175],  cds_values[0] )
  end

  def test_main_parser2
    assert_equal( "Scaffold881", @g3p.gffh["SPU_025302gn"]["SPU_025302-tr"]["mRNA"][0][0] )
#    assert_equal( "Scaffold881", @g3p.agene( "SPU_025302gn" ).atranscript( "SPU_025302-tr" )["mRNA"][0][0] )
    assert_equal( 59053, @g3p.agene( "SPU_025302gn" ).atranscript("SPU_025302-tr").mrna.start )
    assert_equal( 94582, @g3p.agene( "SPU_025302gn" ).atranscript("SPU_025302-tr").mrna.end )
  end

  def test_main_parser3
    assert_equal( 397093, @g3p.scaffoldh["Scaffold881"].size )
  end

  def test_main_parser4
#    assert_equal( "exon", @g3p.gffh["SPU_003298gn"]["SPU_003298-tr"]["exon"][1][2] )
    assert_equal( 324991, @g3p.agene( "SPU_003298gn" ).atranscript( "SPU_003298-tr" ).exons[0].start )
    assert_equal( 94546, @g3p.agene( "SPU_025302gn" ).atranscript( "SPU_025302-tr" ).exons[0].start )
    gidarray = ["SPU_003297gn", "SPU_003298gn", "SPU_011238gn", "SPU_013499gn", "SPU_015325gn", "SPU_017568gn", "SPU_018425gn", "SPU_022816gn", "SPU_022817gn", "SPU_025302gn", "SPU_025908gn", "SPU_025909gn", "SPU_025910gn"]
    assert_equal( gidarray, @g3p.geneids )
  end

=begin
=end

end

