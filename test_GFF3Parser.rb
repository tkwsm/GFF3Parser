#!/usr/bin/ruby

# require 'minitest/unit'
require 'minitest/autorun'
require './GFF3Parser.rb'

class TC_GFF3Parser < MiniTest::Unit::TestCase

  class Test_GFF3Parser
    include GFF3Parser
  end

  def setup
#    @speAgid        = "SPU_025302"
#    @path2speApepfa = "sample_data/spur_protein.fas"
    @path2speAscafa = "sample_data/spur_scaffold.fas"
#    @path2speBscafa = "sample_data/hpul_scaffold.fas"
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
    @g3p = Test_GFF3Parser.new
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
    h     = {}
    ftype = "gene"
    line  = "Scaffold881     GLEAN3-UTR      gene    59053   94582   .       -       .       ID=SPU_025302gn;Name=SPU_025302;"
    assert_equal( ["SPU_025302gn"], @g3p.parse_gff3( line, ftype, h ).keys )
    assert_equal( "Scaffold881",    @g3p.parse_gff3( line, ftype, h )["SPU_025302gn"]["gene"][0] )
    assert_equal( "59053", @g3p.parse_gff3( line, ftype, h )["SPU_025302gn"]["gene"][3] )
    assert_equal( "94582", @g3p.parse_gff3( line, ftype, h )["SPU_025302gn"]["gene"][4] )
    assert_equal( "-",     @g3p.parse_gff3( line, ftype, h )["SPU_025302gn"]["gene"][6] )
  end

  def test_create_scaffold_hash
    assert_equal( 397093, @g3p.create_scaffold_hash( @path2speAscafa )["Scaffold881"].size )
  end

  def test_create_gff_hash
    assert_equal( ["gene", "SPU_025302-tr"], @g3p.create_gff_hash( @path2speAgff )["SPU_025302gn"].keys )
  end

  def test_parse_exon_sites
    exon_values = @g3p.create_gff_hash( @path2speAgff )["SPU_025302gn"]["SPU_025302-tr"]["exon"] 
    assert_equal( ["61717", "62175"], @g3p.parse_exon_sites( exon_values )[0] )
  end

  def test_parse_cds_sites
    cds_values  = @g3p.create_gff_hash( @path2speAgff )["SPU_025302gn"]["SPU_025302-tr"]["exon"]
    assert_equal( ["61717", "62175"], @g3p.parse_cds_sites( cds_values )[0] )
  end

  def test_main_parser2
    assert_equal( "Scaffold881", @g3p.main_parser2( @path2speAgff )["SPU_025302gn"]["SPU_025302-tr"]["mRNA"][0][0] )
  end

  def test_main_parser3
    assert_equal( 397093, @g3p.main_parser3( @path2speAscafa )["Scaffold881"].size )
  end

  def test_main_parser4
    assert_equal( "exon", @g3p.main_parser4( @path2speAgff, @path2speAscafa )["SPU_003298gn"]["SPU_003298-tr"]["exon"][1][2] )
  end


end

