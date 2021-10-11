#!/usr/bin/ruby

require 'rubygems'
# require 'mysql'
require 'bio'
include Bio


module GFF3Parser

  class GFF3Exon

    def initialize( exon_gff3_array, i )
      @exonnumber = (i+1)
      @parenttranscriptid = exon_gff3_array[8]
      @start      = exon_gff3_array[3].to_i
      @end        = exon_gff3_array[4].to_i
      @strand     = exon_gff3_array[6]
    end

    attr_accessor :exonnumber, :start, :end, :strand, 
                  :parenttranscriptid

  end

  class GFF3Intron

    def initialize( intron_gff3_array, i )
      @intronnumber = (i+1)
      @parenttranscriptid = intron_gff3_array[8]
      @start      = intron_gff3_array[3].to_i
      @end        = intron_gff3_array[4].to_i
      @strand     = intron_gff3_array[6]
    end

    attr_accessor :intronnumber, :start, :end, :strand, 
                  :parenttranscriptid

  end

  class GFF3TS

    def initialize( ts_gff3_array )
      @parenttranscriptid = ts_gff3_array[8]
      @start      = ts_gff3_array[3].to_i
      @end        = ts_gff3_array[4].to_i
      @strand     = ts_gff3_array[6]
    end

    attr_accessor :start, :end, :strand, :parenttranscriptid

  end

  class GFF3TE

    def initialize( te_gff3_array )
      @parenttranscriptid = te_gff3_array[8]
      @start      = te_gff3_array[3].to_i
      @end        = te_gff3_array[4].to_i
      @strand     = te_gff3_array[6]
    end

    attr_accessor :start, :end, :strand, :parenttranscriptid

  end

  class GFF3StartCodon

    def initialize( sc_gff3_array )
      @parenttranscriptid = sc_gff3_array[8]
      @start      = sc_gff3_array[3].to_i
      @end        = sc_gff3_array[4].to_i
      @strand     = sc_gff3_array[6]
    end

    attr_accessor :start, :end, :strand, :parenttranscriptid

  end

  class GFF3StopCodon

    def initialize( sc_gff3_array )
      @parenttranscriptid = sc_gff3_array[8]
      @start      = sc_gff3_array[3].to_i
      @end        = sc_gff3_array[4].to_i
      @strand     = sc_gff3_array[6]
    end

    attr_accessor :start, :end, :strand, :parenttranscriptid

  end

  class GFF3CDS

    def initialize( cds_gff3_array, i )
      @cdsnumber = (i+1)
      @parenttranscriptid = cds_gff3_array[8]
      @start      = cds_gff3_array[3].to_i
      @end        = cds_gff3_array[4].to_i
      @strand     = cds_gff3_array[6]
    end

    attr_accessor :cdsnumber, :start, :end, :strand, 
                  :parenttranscriptid

  end

  class GFF3UTR5

    def initialize( utr5_gff3_array, i )
      @exonnumber = (i+1)
      @parenttranscriptid = utr5_gff3_array[8]
      @start      = utr5_gff3_array[3].to_i
      @end        = utr5_gff3_array[4].to_i
      @strand     = utr5_gff3_array[6]
    end

    attr_accessor :utr5number, :start, :end, :strand, 
                  :parenttranscriptid

  end

  class GFF3UTR3

    def initialize( utr3_gff3_array, i )
      @exonnumber = (i+1)
      @parenttranscriptid = utr3_gff3_array[8]
      @start      = utr3_gff3_array[3].to_i
      @end        = utr3_gff3_array[4].to_i
      @strand     = utr3_gff3_array[6]
    end

    attr_accessor :utr3number, :start, :end, :strand, 
                  :parenttranscriptid

  end

  class GFF3mRNA

    def initialize( mrna_start, mrna_end ) 
      @start      = mrna_start
      @end        = mrna_end
    end

    attr_accessor :start, :end

  end

  class GFF3Transcript

    def initialize( transcriptid )

      @transcriptid  = transcriptid
      @parentgeneid  = ""
      @exons         = []
      @exonposisions = []
      @cdss          = []
      @utr5s         = []
      @utr3s         = []
      @introns       = []
      @mrna          = []

    end

    attr_accessor :transcriptid, :parentgeneid, :mrna, :exons, :cdss, 
                  :utr5s, :utr3s, :start_codon, :stop_codon, :introns, 
                  :transcription_start_site, :transcription_end_site

    def sort_exon_sites

      @exons.sort!{|ex, ey| ex.start <=> ey.start }

    end

    def sort_cds_sites

      @cdss.sort!{|cx, cy| cx.start <=> cy.start }

    end

    def parsed_exon_sites
      @exons.collect{|ex| [ex.start, ex.end]}
    end

    def parsed_cds_sites
      @cdss.collect{|cx| [cx.start, cx.end]}
    end

  end

  class GFF3Gene

    def initialize( geneid )

      @geneid          = geneid
      @transcripts     = []
      @transcriptshash = {}

    end

    attr_accessor :geneid, :transcripts

    def atranscript( tid )
      @transcriptshash[ tid ] 
    end

    def set_atranscript( tid, transcript )
      @transcriptshash[ tid ] = transcript
    end

  end

  class GFF3Parsed
  
    include GFF3Parser

    def initialize( in_gff, in_scaffold )
      @exon_sites      = []
      @cds_sites       = []
      @concat_exon_seq = "" 
      @concat_cds_seq  = "" 
      @transcript_seq  = "" 
      @prot_seq        = "" 
      @gffh            = {}
      @scaffoldh       = {}
      @genes           = []
      @geneshash       = {}
      @transcriptshash = {}
      @geneids         = []
      @transcriptids   = []
      @type_array = define_type_array
      main_parser2( in_gff, @gffh )
      main_parser3( in_scaffold, @scaffoldh )
      create_genes
    end

    attr_reader :concat_exon_seq, 
                :concat_cds_seq, :transcript_seq, :cds_seq, 
                :prot_seq, :exon_sites, :cds_sites

    attr_accessor :gffh, :scaffoldh, :genes, :geneids, :transcriptids

    def create_genes 
      @geneids = @gffh.keys.sort
      @gffh.each_key{ |gid| @genes << GFF3Gene.new( gid ) }
      @genes.each do |gene|
        @geneshash[ gene.geneid ] = gene
        @gffh[ gene.geneid ].keys.each do |tid|
          next if tid == "gene"
          transcript = GFF3Transcript.new( tid )
          gene.set_atranscript( tid, transcript )
          transcript.parentgeneid = gene.geneid
          gid = gene.geneid
          @type_array.each do |ftype|
            data_array = @gffh[gid][tid][ftype]
            if    ftype == "mRNA"
              mrna_start = data_array[0][3].to_i
              mrna_end   = data_array[0][4].to_i
              transcript.mrna         = GFF3mRNA.new( mrna_start, mrna_end )
            elsif ftype == "exon"
              data_array.each_with_index do |aexon, i|
                transcript.exons << GFF3Exon.new( aexon, i ) 
                transcript.cdss << GFF3CDS.new( aexon, i )   ###
              end
            elsif ftype == "CDS"
              data_array.each_with_index do |acds, i|
                transcript.cdss << GFF3CDS.new( acds, i )
              end
            elsif ftype == "5-utr" or ftype == "five_prime_UTR"
              data_array.each_with_index do |a5utr, i|
                transcript.utr5s << GFF3UTR5.new( a5utr, i )
              end
            elsif ftype == "3-utr" or ftype == "three_prime_UTR"
              data_array.each_with_index do |a3utr, i|
                transcript.utr3s << GFF3UTR3.new( a3utr, i )
              end
            elsif ftype == "intron"
              data_array.each_with_index do |intron, i|
                transcript.introns << GFF3Intron.new( intron, i )
              end
            elsif ftype == "start_codon"
              transcript.start_codon = GFF3StartCodon.new( data_array )
            elsif ftype == "stop_codon"
              transcript.stop_codon  = GFF3StopCodon.new( data_array )
            elsif ftype == "transcription_start_site"
              transcript.transcription_start_site = GFF3TS.new( data_array )
            elsif ftype == "transcription_end_site"
              transcript.transcription_end_site   = GFF3TE.new( data_array )
            else
              p [ tid, ftype ]; p "something strange! 1"; exit;
            end
          end
          gene.transcripts << transcript
        end
      end
    end

    def agene( gid )
      @geneshash[ gid ] 
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
        STDERR.puts "Create trancript_seq after direction was initialized"
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

  end # Class GFF3Parsed


  def sampletest
    "good result"
  end

  def print_usage
    "\n ./GFF3Parser.rb --gff      <gff3> \n \
                        --scaffold <project_id> \n\n "
  end

  def create_gene_hash

    h = {}
    type_array = [ "gene" ]
    type_array.each do |type|
      h[type] = []
    end

    return h

  end # create_gene_hash

  def define_type_array 
    type_array = %w( CDS exon intron start_codon stop_codon mRNA transcription_end_site transcription_start_site 3-utr 5-utr)
    return type_array 
  end

  def create_transcript_hash
  
    h = {}
    type_array = define_type_array
    type_array.each do |type|
      h[type] = []
    end

    return h

  end # create_transcript_hash

  def parse_gff3( line, ftype, h )

    dataline   = []
    seqid      = "" ;  
    source     = ""; 
    sstart     = 0; 
    send       = 0; 
    score      = ".";
    strand     = "+"; 
    phase      = "."; 
    attributes = "";
    gid        = ""; 
    tid        = "";
    type_array = define_type_array
  
    a = line.chomp.split("\s")
    seqid  = a[0]
    source = a[1]
    ftype  = ftype
    sstart = a[3]
    send   = a[4]
    score  = a[5]
    strand = a[6]
    phase  = a[7]
    attributes = a[8..-1].join(" ")
    dataline = [seqid, source, ftype, sstart, send, score, strand, phase ]
  
    if ftype == "gene" and attributes =~ /ID/
  
      gid = attributes.slice(/ID=([^\;]+)/, 1)
      if gid == nil
        STDERR.puts "gid is strange 5"
        STDERR.puts attributes
      end
  
      h[gid] = create_gene_hash if h[gid] == nil
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

      h[gid]      = create_gene_hash       if h[gid] == nil
      h[gid][tid] = create_transcript_hash if h[gid][tid] == nil
      h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]
  
    elsif ftype == "exon" and attributes =~ /Parent/

      tids = attributes.slice(/Parent=([^\;]+)/, 1)
      if    tids =~ /\,/
        tids.split(",").each do |tid|
          gid = tid.slice(/=?(\S+)-mRNA/, 1) 
          h[gid]      = create_gene_hash       if h[gid] == nil
          h[gid][tid] = create_transcript_hash if h[gid][tid] == nil
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

        h[gid]      = create_gene_hash       if h[gid] == nil
        h[gid][tid] = create_transcript_hash if h[gid][tid] == nil
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

      h[gid]      = create_gene_hash       if h[gid] == nil
      h[gid][tid] = create_transcript_hash if h[gid][tid] == nil
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

      h[gid]      = create_gene_hash       if h[gid] == nil
      h[gid][tid] = create_transcript_hash if h[gid][tid] == nil
      h[gid][tid][ ftype ] << dataline + ["ID=#{tid}\;Parent=#{gid}"]
  
    else
      p "something bad : #{ftype} : #{attributes} : #{line} : #{a}4 "
      exit
  
    end
  
    return h 
  
  end # parse_gff3( line, ftype, h )

  def create_scaffold_hash( input_scaffold, scaffold_h )
    
    ff = FlatFile.new(FastaFormat, open(input_scaffold) )
    ff.each do |e|
      scaffold_id = e.definition.slice(/(\S+)/, 1)
      scaffold_h[ scaffold_id ] = e.naseq
    end
    return scaffold_h

  end

  def create_gff_hash( input_gff, h )

    scaffold_id = ""
    transcript_id = ""
    cdna_seq = ""
    prot_seq = ""
    cdna = ""
    a = []
    ftype = ""
    open( input_gff ).each_with_index do |x, i|
      a = x.chomp.split("\t")
      ftype = a[2]
      next if x =~ /^\#/ 
      next if x == "\n"

      if ftype == "gene"  
        h = parse_gff3( x, "gene", h )
      elsif ftype == "transcript" or ftype ==  "mRNA"
        h = parse_gff3( x, "mRNA", h )
      elsif ftype == "exon"  
        h = parse_gff3( x, "exon", h )
      elsif ftype == "CDS"  
        h = parse_gff3( x, "CDS", h )
      elsif ftype == "start_codon"
        h = parse_gff3( x, "start_codon", h )
      elsif ftype == "stop_codon"
        h = parse_gff3( x, "stop_codon", h )
      elsif ftype == "3-utr"
        h = parse_gff3( x, "3-utr", h )
      elsif ftype == "three_prime_UTR"
        h = parse_gff3( x, "3-utr", h )
      elsif ftype == "5-utr"
        h = parse_gff3( x, "5-utr", h )
      elsif ftype == "five_prime_UTR"
        h = parse_gff3( x, "5-utr", h )
      elsif ftype == "transcription_start_site"
        h = parse_gff3( x, "transcription_start_site", h )
      elsif ftype == "transcription_end_site"
        h = parse_gff3( x, "transcription_end_site", h )
      elsif ftype == "intron"
      else
        p x
        p "something strange! 1"
        exit
      end
    end  # open( input_gff )
  end

  def parse_exon_sites( exon_values )
    exon_sites = []
    exon_values.each do |e|
      exon_sites << e[3..4]
    end
    exon_sites.sort!{|ex, ey| ex[0].to_i <=> ey[0].to_i}
    return exon_sites
  end

  def parse_cds_sites( cds_values )

    if    cds_values[0].strand == "+"
      cds_values.sort!{|ex, ey| ex.start <=> ey.start }
    elsif cds_values[0].strand == "-"
      cds_values.sort!{|ex, ey| ey.end   <=> ex.end   }
    end
    return cds_values
  end

  def main_parser2( in_gff, h)

    h = create_gff_hash( in_gff, h )
    return h

  end

  def main_parser3( in_scaffold, scah )

    scah = create_scaffold_hash( in_scaffold, scah )
    return scah

  end

  def main_parser( in_gff, in_scaffold, out_cdna, out_prot, out_cdsnucl )

    scah = GFF3Parser.create_scaffold_hash( in_scaffold )

    h.each_key do |gid|
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

