use strict;
use warnings;
    
## FUNCOES

sub GetORF{ # esta função chama um programa do pacote EMBOSS para encontrar ORFS a partir de um arquivo fasta
    my $INFILE = $_[0];
    # o parametro "-find 3" define que as ORFs encontradas são retornadas como nucleotideos
    `getorf -sequence $INFILE -outseq $INFILE.orf.temp -find 3`
}

sub FilterORF{
    my $INFILE = $_[0];
    my $min_length = 1;
    `perl src/FilterFastaLength.pl $min_length $INFILE > $INFILE.filtered.temp`;
    `cp $INFILE.filtered.temp $INFILE`;
}

sub DiFreq{ # Função para calcular a frequencia de di e tri nucleotideos, e normalizando pela razão do esperado para cada
    my $ORF_fasta = $_[0];
    my ($FASTA, $secondlast, $last);
    my ($len, $i, $nt, $nt_pair, $nt_trio);
    my (%mono_nt, %di_nt, %tri_nt);
    my ($each_mono_nt, $each_di_nt ,$each_tri_nt);
    my ($first_nt, $second_nt, $third_nt);
    my ($obs_di, $expec_di, $obs_tri, $expec_tri);
    my (%odds_ratio_dint, %odds_ratio_trint);
    my @Seq = ();
    my $temp_seq = "";
    # Abrindo o arquivo no modo de leitura no ponteiro "$FASTA"
    open $FASTA, '<', $ORF_fasta or die $!; # $! imprime o erro apropriado
    while (my $line = <$FASTA>) {
        chomp($line);
        if ($line =~ /^>+/){ 
            if ($temp_seq ne ""){ # a sequencia anterior é adicionada quando o proximo ">" é encontrado
                push @Seq, $temp_seq;
                $temp_seq = ""; # apaga as linhas acumuladas
            }
            next; # Para sair do laço e não atualizar $temp_seq
        }
        $temp_seq .= $line; #concatena o conteudo de $line, caso não tenha encontrado ">"
    }
    if ($temp_seq ne ""){ # adiciona a sequencia anterior no array @Seq
        push @Seq, $temp_seq;
        $temp_seq = "";
    }
    my $total_len = 0;
    foreach my $seq (@Seq){
        $len = length $seq; # salva o tamanho da sequencia
        $total_len = $total_len + $len;
        for $i (0..$len-3) {  # laço percorrendo a sequencia até o posição -3 por causa da contagem de trint
            $each_mono_nt = substr($seq, $i, 1);
            $each_di_nt = substr($seq, $i, 2);
            $each_tri_nt = substr($seq, $i, 3);
            #print $each_mono_nt," ",$each_di_nt, " ", $each_tri_nt,"\n";
            $mono_nt{$each_mono_nt}++;
            $di_nt{$each_di_nt}++;
            $tri_nt{$each_tri_nt}++;
            $last = uc substr($seq, -1); # ultima letra da seq
            $secondlast = uc substr($seq, -2,-1); #penultima letra
        }
        $mono_nt{$secondlast}++; # incrementa o ultimo dinucleotideo, que nao foi comtemplado no laço anterior
        $mono_nt{$last}++; # incrementa o ultimo nt
        $di_nt{"$secondlast$last"}++; # incrementa o ultimo dint
    }
    #print "-"x30, "\nSingle nucleotide frequency:\n";
    for $nt (sort keys %mono_nt) {
        #print "$nt\t", $mono_nt{$nt} / $total_len, "\n";
    }
    #print "\n", "-"x30, "\nDinucleotide frequency:\nDinucleotide\tObs. freq.\tExp. freq.\tOdds ratio:\n";
    for $nt_pair (sort keys %di_nt) {
        ($first_nt, $second_nt) = split //, $nt_pair; # Divide a string da key em dois caracteres
        $obs_di = $di_nt{$nt_pair} / ($total_len-1);
        $expec_di = $mono_nt{$first_nt} * $mono_nt{$second_nt} / $total_len / $total_len;
        $odds_ratio_dint{"$nt_pair"} = ($obs_di / $expec_di);
        #print "$nt_pair\t", $obs_di, "\t", $expec_di, "\t", $odds_ratio_dint{$nt_pair}, "\n"; # valor esperado, baseado na frequencia dos single nt
    }
    return %odds_ratio_dint; 
}

sub TriFreq{ # Função para calcular a frequencia de di e tri nucleotideos, e normalizando pela razão do esperado para cada
    my $ORF_fasta = $_[0];
    my ($FASTA, $secondlast, $last);
    my ($len, $i, $nt, $nt_pair, $nt_trio);
    my (%mono_nt, %di_nt, %tri_nt);
    my ($each_mono_nt, $each_di_nt ,$each_tri_nt);
    my ($first_nt, $second_nt, $third_nt);
    my ($obs_di, $expec_di, $obs_tri, $expec_tri);
    my (%odds_ratio_dint, %odds_ratio_trint);
    my @Seq = ();
    my $temp_seq = "";
    # Abrindo o arquivo no modo de leitura no ponteiro "$FASTA"
    open $FASTA, '<', $ORF_fasta or die $!; # $! imprime o erro apropriado
    while (my $line = <$FASTA>) {
        chomp($line);
        if ($line =~ /^>+/){ 
            if ($temp_seq ne ""){ # a sequencia anterior é adicionada quando o proximo ">" é encontrado
                push @Seq, $temp_seq;
                $temp_seq = ""; # apaga as linhas acumuladas
            }
            next; # Para sair do laço e não atualizar $temp_seq
        }
        $temp_seq .= $line; #concatena o conteudo de $line, caso não tenha encontrado ">"
    }
    if ($temp_seq ne ""){ # adiciona a sequencia anterior no array @Seq
        push @Seq, $temp_seq;
        $temp_seq = "";
    }
    my $total_len = 0;
    foreach my $seq (@Seq){
        $len = length $seq; # salva o tamanho da sequencia
        $total_len = $total_len + $len;
        for $i (0..$len-3) {  # laço percorrendo a sequencia até o posição -3 por causa da contagem de trint
            $each_mono_nt = substr($seq, $i, 1);
            $each_di_nt = substr($seq, $i, 2);
            $each_tri_nt = substr($seq, $i, 3);
            #print $each_mono_nt," ",$each_di_nt, " ", $each_tri_nt,"\n";
            $mono_nt{$each_mono_nt}++;
            $di_nt{$each_di_nt}++;
            $tri_nt{$each_tri_nt}++;
            $last = uc substr($seq, -1); # ultima letra da seq
            $secondlast = uc substr($seq, -2,-1); #penultima letra
        }
        $mono_nt{$secondlast}++; # incrementa o ultimo dinucleotideo, que nao foi comtemplado no laço anterior
        $mono_nt{$last}++; # incrementa o ultimo nt
        $di_nt{"$secondlast$last"}++; # incrementa o ultimo dint
    }
    #print "-"x30, "\nSingle nucleotide frequency:\n";
    for $nt (sort keys %mono_nt) {
        #print "$nt\t", $mono_nt{$nt} / $total_len, "\n";
    }
    #print "\n", "-"x30, "\nTrinucleotide frequency:\nTrinucleotide\tObs. freq.\tExp. freq.\tOdds ratio:\n";
    for $nt_trio (sort keys %tri_nt){
        ($first_nt, $second_nt, $third_nt) = split //, $nt_trio; # Divide a string da key em tres caracteres
        $obs_tri = $tri_nt{$nt_trio} / ($total_len-2);
        $expec_tri = $mono_nt{$first_nt} * $mono_nt{$second_nt} * $mono_nt{$third_nt} / $total_len / $total_len / $total_len;
        $odds_ratio_trint{"$nt_trio"} = ($obs_tri / $expec_tri);
        #print "$nt_trio\t", $obs_tri, "\t", $expec_tri, "\t", $odds_ratio_trint{$nt_trio}, "\n"; # valor esperado, baseado na frequencia dos single nt
    }
    return %odds_ratio_trint; 
}

sub SaveTabFiles{
    my ($INFILE, $reference_di, $reference_tri) = @_;
    my %odds_ratio_dint = %{ $reference_di }; #derreferenciando
    my %odds_ratio_trint = %{ $reference_tri }; #derreferenciando
    my $OUTFILE_di = $INFILE."DiFreq.temp";
    my $OUTFILE_tri = $INFILE."TriFreq.temp";
    my $OUTFILE_concat = $INFILE."ConcatFreq.temp";
    open (my $fh_di, '>', $OUTFILE_di) or die "ERROR opening $OUTFILE_di\n";
    open (my $fh_tri, '>', $OUTFILE_tri) or die "ERROR opening $OUTFILE_tri\n";
    open (my $fh_concat, '>', $OUTFILE_concat) or die "ERROR opening $OUTFILE_concat\n";
    for my $nt_pair (sort keys %odds_ratio_dint) {
        (my $first_nt, my $second_nt) = split //, $nt_pair; # Divide a string da key em dois caracteres
        print $fh_di "$nt_pair\t", $odds_ratio_dint{$nt_pair}, "\n";
        print $fh_concat "$nt_pair\t", $odds_ratio_dint{$nt_pair}, "\n";
    }
    for my $nt_trio (sort keys %odds_ratio_trint){
        (my $first_nt, my $second_nt, my $third_nt) = split //, $nt_trio; # Divide a string da key em tres caracteres
        print $fh_tri "$nt_trio\t", $odds_ratio_trint{$nt_trio}, "\n";
        print $fh_concat "$nt_trio\t", $odds_ratio_trint{$nt_trio}, "\n";
    }
    close $fh_di; 
    close $fh_tri;
    close $fh_concat;
    print "Perfis de di- e tri- nts salvos com sucesso.\n";
}

sub CreateDatabase{
    my $db_file = $_[0];
    my @nt = ('A', 'C', 'G', 'T');
    print "Criando banco de dados.\n";
    open(my $fh_db, '>', $db_file);
    print $fh_db "lib_name\n";
    foreach my $nt1 (@nt){
        foreach my $nt2 (@nt){
            print $fh_db $nt1, $nt2, "\n";
        }
    }
    foreach my $nt1 (@nt){
        foreach my $nt2 (@nt){
            foreach my $nt3 (@nt){
                print $fh_db $nt1, $nt2, $nt3, "\n";
            }
        }
    }
    print "Banco de dados $db_file criado com sucesso.\n";
}

sub AddDatabase{
    my $db_name = $_[0];
    my $lib_name = $_[1];
    my $reference_di = $_[2];
    my $reference_tri = $_[3];
    my %odds_ratio_dint = %{ $reference_di }; #derreferenciando
    my %odds_ratio_trint = %{ $reference_tri }; #derreferenciando
   
    my $db_file = "data/$db_name.db";
    my @nt = ('A', 'C', 'G', 'T');
    if (-f $db_file){ ## testa se o arquivo já existe, caso não exista, formata os campos do banco de dados 
        print "Banco de dados: $db_file existente.\n";
    }
    else{
        print $db_file,"\n";
        CreateDatabase( $db_file );
    }
    open(my $fh_data, '>', "data/$lib_name.-alldata.scaled.temp");
    print $fh_data $lib_name, "\n";
    foreach my $nt1 (@nt){
        foreach my $nt2 (@nt){
            if (exists $odds_ratio_dint{"$nt1$nt2"}){
                print $fh_data $odds_ratio_dint{"$nt1$nt2"}, "\n";
            }
            else {
                print $fh_data "0\n";
            }
        }
    }
    foreach my $nt1 (@nt){
        foreach my $nt2 (@nt){
            foreach my $nt3 (@nt){
                if (exists $odds_ratio_trint{"$nt1$nt2$nt3"}){
                    print $fh_data $odds_ratio_trint{"$nt1$nt2$nt3"}, "\n";
                }
                else{
                    print $fh_data "0\n"
                }
            }
        }
    }
    close($fh_data);
    `paste -d'\t' $db_file data/$lib_name.-alldata.scaled.temp > $db_file.temp`;
    `cp $db_file.temp $db_file`;
    print "Banco de dados atualizado com sucesso!\n";
}

sub PlotProfile{
    my $db_file = $_[0];
    my $lib = $_[1];
    `Rscript --no-save src/PlotProfile.R $db_file $lib &`;
}

sub PearsonCorrelation{ 
    my $db_file = $_[0];
    `Rscript --no-save src/PearsonCorrelation.R $db_file &`;
}

return 1;  #retorna TRUE quando executado corretamente
