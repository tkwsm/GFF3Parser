#!/usr/bin/ruby

require 'rubygems'
# require 'mysql'
require 'bio'
include Bio

class GFF3Parsed

  def initialize
    @gid = ""
    @tid = ""
    @scaffold_id     = ""
    @direction       = ""
    @exon_sites      = []
    @cds_sites       = []
    @concat_exon_seq = "" 
    @concat_cds_seq  = "" 
    @transcript_seq  = "" 
    @prot_seq        = "" 
  end

  attr_reader :gid, :tid, :scaffold_id, :direction, :concat_exon_seq, :concat_cds_seq, :transcript_seq, :cds_seq, :prot_seq, :exon_sites, :cds_sites

  def get_gid( gid )
    @gid = gid
  end

  def get_tid( tid )
    @tid = tid
  end

  def get_scaffold_id( scaffold_id )
    @scaffold_id = scaffold_id
  end

  def get_direction( direction )
    @direction = direction
  end

  def get_exon_sites( exon_sites )
    exon_sites.each{|e| @exon_sites << [e[0].to_i, e[1].to_i]}
    @exon_sites.sort!{|x, y| x[0].to_i <=> y[0].to_i}
  end

  def get_cds_sites( cds_sites )
    cds_sites.each{|e| @cds_sites << [e[0].to_i, e[1].to_i]}
    @cds_sites.sort!{|x, y| x[0].to_i <=> y[0].to_i}
  end


  def parse_nucl_sequence( dna_seq )
    concat_exon_seq = ""
    formatted_dna_seq = Sequence::NA.new( dna_seq )
    @exon_sites.each do |exon|
      exon_start, exon_end = exon
      concat_exon_seq += formatted_dna_seq.subseq( exon_start, exon_end )
    end
    @concat_exon_seq = concat_exon_seq
    if    @direction == "+"
      @transcript_seq = Sequence::NA.new( concat_exon_seq )
    elsif @direction == "-"
      @transcript_seq = Sequence::NA.new( concat_exon_seq ).complement
    else
      @transcript_seq = Sequence::NA.new( concat_exon_seq )
      STDERR.puts "please create trancript_seq after direction was initialized"
    end
  end

  def parse_cds_sequence( dna_seq )
    concat_cds_seq = ""
    formatted_dna_seq = Sequence::NA.new( dna_seq )
    @cds_sites.each do |cds|
      cds_start, cds_end = cds
      concat_cds_seq += formatted_dna_seq.subseq( cds_start, cds_end )
    end
    @concat_cds_seq = concat_cds_seq
    if    @direction == "+"
      @cds_seq = Sequence::NA.new( concat_cds_seq )
    elsif @direction == "-"
      @cds_seq = Sequence::NA.new( concat_cds_seq ).complement
    else
      @cds_seq = Sequence::NA.new( concat_cds_seq )
      STDERR.puts "please create cds_seq after direction was initialized"
    end
    @protein_seq = @cds_seq.translate
  end

  def show_transcript_sequence
    @transcript_seq.to_fasta( @tid, 60 )
  end

  def show_protein_sequence
    @protein_seq.to_fasta( @tid, 60 )
  end

  def show_cds_sequence
    @cds_seq.to_fasta( @tid, 60 )
  end

end

module GFF3Parser

  def sampletest
    "good result"
  end

  def print_usage
    "\n ./GFF3Parser.rb --gff      <gff3> \n \
                        --scaffold <project_id> \n\n "
  end

  def GFF3Parser.create_gene_hash

    h = {}
    type_array = [ "gene" ]
    type_array.each do |type|
      h[type] = []
    end

    return h

  end # create_gene_hash


  def GFF3Parser.create_transcript_hash
  
    h = {}
    type_array = %w( CDS exon intron start_codon stop_codon mRNA transcription_end_site transcription_start_site 3-utr 5-utr)
    type_array.each do |type|
      h[type] = []
    end

    return h

  end # create_transcript_hash


  def GFF3Parser.sort_features( features_array )
  
    if features_array == []
    elsif features_array[0][6] == "+"
      features_array.sort!{|x, y| x[3].to_i <=> y[3].to_i }
    elsif features_array[0][6] == "-"
      features_array.sort!{|x, y| y[3].to_i <=> x[3].to_i }
    else
      p "something wrong #{features_array[0]} 2"
      p "something wrong #{features_array[0][6]} 3"
      exit
    end
    return features_array

  end # sort_features( features_array )

  def GFF3Parser.parse_gff3( line, ftype, h )

    dataline = []
    seqid  = "" ;  
    source = ""; 
    sstart = 0; 
    send   = 0; 
    score = ".";
    strand = "+"; 
    phase  = "."; 
    attributes = "";
    gid    = ""; 
    tid = "";
    type_array = %w( CDS exon intron start_codon stop_codon transcription_end_site transcription_start_site 5-utr 3-utr )
  
    a = line.chomp.split("\t")
    seqid  =  a[0]
    source = a[1]
    ftype  = ftype
    sstart = a[3]
    send   = a[4]
    score  = a[5]
    strand = a[6]
    phase  = a[7]
    attributes = a[8]
    dataline = [seqid, source, ftype, sstart, send, score, strand, phase ]
  
    if ftype == "gene" and attributes =~ /ID/
  
      gid = attributes.slice(/ID=([^\;]+)/, 1)
      if gid == nil
        STDERR.puts "gid is strange 5"
        STDERR.puts attributes
      end
  
      h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
      h[gid][ ftype ] = dataline + ["ID=#{gid}\;Name=#{gid}"]
  
    elsif ftype == "mRNA" and attributes =~ /Parent/
  
      if attributes =~ /Parent=[^\;]+gn/
        gid = attributes.slice(/Parent=([^\;]+gn)/, 1)
      else
        gid = attributes.slice(/Parent=([^\;]+)/, 1)
      end
      tid = attributes.slice(/ID=([^\;]+)/, 1)
      if gid == nil
        STDERR.puts "gid is strange 6"
        STDERR.puts attributes
      end

      h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
      h[gid][tid] = GFF3Parser.create_transcript_hash if h[gid][tid] == nil
      h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]
  
    elsif ftype == "exon" and attributes =~ /Parent/

      tids = attributes.slice(/Parent=([^\;]+)/, 1)
      if    tids =~ /\,/
        tids.split(",").each do |tid|
          gid = tid.slice(/=?(\S+)-mRNA/, 1) 
          h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
          h[gid][tid] = GFF3Parser.create_transcript_hash if h[gid][tid] == nil
          h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]
        end
      else
        if tids =~ /-mRNA/
          tid = tids
          gid = tid.slice(/=?(\S+)-mRNA/, 1) 
        elsif tids =~ /\S\.\S-tr/
          tid = tids
          gid = tid.slice(/=?([^\.]+)\.\S-tr/, 1)  + "gn"
        elsif tids =~ /\S-tr/
          tid = tids
          gid = tid.slice(/=?([^\.]+)-tr/, 1)  + "gn"
        end

        h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
        h[gid][tid] = GFF3Parser.create_transcript_hash if h[gid][tid] == nil
        h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]
      end

    elsif ftype == "CDS" and attributes =~ /Parent/
  
      tid = attributes.slice(/Parent=([^\;]+)/, 1)
      if tid =~ /-mRNA/
        gid = tid.slice(/=?(\S+)-mRNA/, 1) 
      elsif tid =~ /\S\.\S-tr/
        gid = tid.slice(/=?([^\.]+)\.\S-tr/, 1)  + "gn"
      elsif tid =~ /\S-tr/
        gid = tid.slice(/=?([^\.]+)-tr/, 1)  + "gn"
      else
        gid = tid
      end

      h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
      h[gid][tid] = GFF3Parser.create_transcript_hash if h[gid][tid] == nil
      h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]

    elsif type_array.include?( ftype ) and attributes =~ /Parent/
  
      tid = attributes.slice(/Parent=([^\;]+)/, 1)
      if tid =~ /-mRNA/
        gid = tid.slice(/=?(\S+)-mRNA/, 1) 
      elsif tid =~ /\S\.\S-tr/
        gid = tid.slice(/=?([^\.]+)\.\S-tr/, 1)  + "gn"
      elsif tid =~ /\S-tr/
        gid = tid.slice(/=?([^\.]+)-tr/, 1)  + "gn"
      else
        gid = tid
      end

      h[gid] = GFF3Parser.create_gene_hash if h[gid] == nil
      h[gid][tid] = GFF3Parser.create_transcript_hash if h[gid][tid] == nil
      h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]
  
    else
      p "something bad #{line} 4"
      exit
  
    end
  
    return h 
  
  end # parse_gff3( line, ftype, h )

  def GFF3Parser.create_scaffold_hash( input_scaffold )
    
    scaffold_h = {}
    ff = FlatFile.new(FastaFormat, open(input_scaffold) )
    ff.each do |e|
      scaffold_id = e.definition.slice(/(\S+)/, 1)
      scaffold_h[ scaffold_id ] = e.naseq
    end
    return scaffold_h

  end

  def GFF3Parser.create_gff_hash( input_gff )

    scaffold_id = ""
    transcript_id = ""
    cdna_seq = ""
    prot_seq = ""
    cdna = ""
    a = []
    h = {}
    ftype = ""
    open( input_gff ).each_with_index do |x, i|
      a = x.chomp.split("\t")
      ftype = a[2]
      next if x =~ /^\#/ 
      next if x == "\n"

      if ftype == "gene"  
        h = GFF3Parser.parse_gff3( x, "gene", h )
      elsif ftype == "transcript" or ftype ==  "mRNA"
        h = GFF3Parser.parse_gff3( x, "mRNA", h )
      elsif ftype == "exon"  
        h = GFF3Parser.parse_gff3( x, "exon", h )
      elsif ftype == "CDS"  
        h = GFF3Parser.parse_gff3( x, "CDS", h )
      elsif ftype == "start_codon"
        h = GFF3Parser.parse_gff3( x, "start_codon", h )
      elsif ftype == "stop_codon"
        h = GFF3Parser.parse_gff3( x, "stop_codon", h )
      elsif ftype == "3-utr"
        h = GFF3Parser.parse_gff3( x, "3-utr", h )
      elsif ftype == "three_prime_UTR"
        h = GFF3Parser.parse_gff3( x, "3-utr", h )
      elsif ftype == "5-utr"
        h = GFF3Parser.parse_gff3( x, "5-utr", h )
      elsif ftype == "five_prime_UTR"
        h = GFF3Parser.parse_gff3( x, "5-utr", h )
      elsif ftype == "transcription_start_site"
        h = GFF3Parser.parse_gff3( x, "transcription_start_site", h )
      elsif ftype == "transcription_end_site"
        h = GFF3Parser.parse_gff3( x, "transcription_end_site", h )
      elsif ftype == "intron"
      else
        p x
        p "something strange! 1"
        exit
      end
    end  # open( input_gff )
    return h
  end

  def GFF3Parser.parse_exon_sites( exon_values )
    exon_sites = []
    exon_values.each do |e|
      exon_sites << e[3..4]
    end
    return exon_sites
  end

  def GFF3Parser.parse_cds_sites( cds_values )
    cds_sites = []
    cds_values.each do |e|
      cds_sites << e[3..4]
    end
    return cds_sites
  end

  def GFF3Parser.main_parser2( in_gff )

    gffh = GFF3Parser.create_gff_hash( in_gff )
    return gffh

  end

  def GFF3Parser.create_scagffh( gffh )

    scagffh = {}
    gffh.each_key do |gid|
      scaid = gffh[gid]["gene"][0]
      scagffh[ scaid ] = [] if scagffh[ scaid ] == nil
      scagffh[ scaid ] << gid
    end
    return scagffh

  end

  def GFF3Parser.main_parser3( in_scaffold )

    scah = GFF3Parser.create_scaffold_hash( in_scaffold )
    return scah

  end

  def GFF3Parser.main_parser4( in_gff, in_scaffold )

    gffh = GFF3Parser.create_gff_hash( in_gff )
    scah = GFF3Parser.create_scaffold_hash( in_scaffold )

    gffh.each_key do |gid|
      gffh[gid].each_key do |tid|
        next if tid == "gene"
        cdnas = GFF3Parsed.new
        cdnas.get_tid( tid )

        exon_values = gffh[gid][tid]["exon"]
        cds_values = gffh[gid][tid]["CDS"]

        scaffold_id = gffh[gid][tid]["mRNA"][0][0]
        direction   = gffh[gid][tid]["mRNA"][0][6]

        exons = GFF3Parser.parse_exon_sites( exon_values ) 
        cdss  = GFF3Parser.parse_cds_sites( cds_values ) 

        cdnas.get_exon_sites( exons )
        cdnas.get_cds_sites( cdss )
        cdnas.get_scaffold_id( scaffold_id )
        cdnas.get_direction( direction )
        dna_seq = scah[ scaffold_id ]
        cdnas.parse_nucl_sequence( dna_seq )
        cdnas.parse_cds_sequence( dna_seq )
      end
    end

  end  #  main_parser2

  def GFF3Parser.main_parser( in_gff, in_scaffold, out_cdna, out_prot, out_cdsnucl )

    gffh = GFF3Parser.create_gff_hash( in_gff )
    scah = GFF3Parser.create_scaffold_hash( in_scaffold )

    gffh.each_key do |gid|
      gffh[gid].each_key do |tid|
        next if tid == "gene"
        cdnas = GFF3Parsed.new
        cdnas.get_tid( tid )

        exon_values = gffh[gid][tid]["exon"]
        cds_values = gffh[gid][tid]["CDS"]

        scaffold_id = gffh[gid][tid]["mRNA"][0][0]
        direction   = gffh[gid][tid]["mRNA"][0][6]

        exons = GFF3Parser.parse_exon_sites( exon_values ) 
        cdss  = GFF3Parser.parse_cds_sites( cds_values ) 

        cdnas.get_exon_sites( exons )
        cdnas.get_cds_sites( cdss )
        cdnas.get_scaffold_id( scaffold_id )
        cdnas.get_direction( direction )
        dna_seq = scah[ scaffold_id ]
        cdnas.parse_nucl_sequence( dna_seq )
        cdnas.parse_cds_sequence( dna_seq )
        out_cdna.print cdnas.show_transcript_sequence
        out_prot.print cdnas.show_protein_sequence
        out_cdsnucl.print cdnas.show_cds_sequence
      end
    end

  end  #  main_parser

end # module GFF3Parser

###########################################################
if $0 == __FILE__
###########################################################

  if ARGV.size == 0
    print "./program <gff> <scaffold>\n"
    print GFF3Parser.print_usage
    exit
  end

  in_gff      = ARGV.shift
  in_scaffold = ARGV.shift
  in_gff_name = in_gff.split("/")[-1]
  out_cdna    = File.new( "#{in_gff_name}.cdna", "w+")
  out_prot    = File.new( "#{in_gff_name}.prot", "w+")
  out_cdsnucl = File.new( "#{in_gff_name}.cds", "w+")

  GFF3Parser.main_parser( in_gff, in_scaffold, out_cdna, out_prot, out_cdsnucl)

###########################################################
end  # __FILE__
###########################################################
#

