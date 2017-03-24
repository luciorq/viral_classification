#!/usr/bin/perl

## Programa criado como parte do projeto final proposto pela disciplina
## Ambientes de Computação - UFMG, realizada no segundo semestre de 2016

## Autores: Lúcio Rezende Queiroz & Elisson Nogueira Lopes

## Este programa encontra a frequencia com que dinucleotideos e
## trinucleotideos ocorrem em todas as sequencias de um dado arquivo FASTA,
## Comparando com as frequencias esperadas e criando um perfil para grupo,
## comparando-os com um agrupamento hierarquico, bsaseado em correlação de Pearson.

## Este programa depende do programa getorf do PACOTE EMBOSS, disponivel em: 
## http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html

## Utilização: perl DiTriNucFreq.pl <sequences.fasta> <database_name>
## Exemplo: 

#bibliotecas
use strict;
use warnings;
use src::functions;

#variaveis declaradas
my ($FASTA, $ORF_fasta);
my (%odds_ratio_dint, %odds_ratio_trint);
my ($db_name, $db_file, @lib_name, $entry_name);

# PROGRAMA PRINCIPAL

print "\n","-"x60,"\n","Executando para: ",$ARGV[0],"\n";

GetORF( $ARGV[0] ); #Encontra todas as ORFs em um arquivo FASTA
$ORF_fasta = $ARGV[0].".orf.temp";
FilterORF( $ORF_fasta ); #Filtra as sequencias encontradas pela função GetORF por tamanho
if (-z $ORF_fasta ) { #Testa se o arquivo está vazio, encerrando a execução para este arquivo caso esteja
	die "Arquivo ",$ORF_fasta," vazio ou sem ORFs válidas.\n";
}

%odds_ratio_dint = DiFreq($ORF_fasta); # calcula frequencia de Dinucleotideos e faz normalização por "odds ratio"
%odds_ratio_trint = TriFreq($ORF_fasta);  # calcula frequencia de Trinucleotideos e faz normalização por "odds ratio"

#print (%odds_ratio_dint,"\n\n" ,%odds_ratio_trint,"\n");

SaveTabFiles( $ORF_fasta, \%odds_ratio_dint, \%odds_ratio_trint ); # salva o arquivo intermediario com os valores de odds ratio

$db_name = $ARGV[1];
$db_file = "data/$db_name.db";
@lib_name = split(/[.]/ ,$ORF_fasta);
@lib_name = split(/[\/]/, $lib_name[0]);
$entry_name = $lib_name[1];

AddDatabase( $db_name, $entry_name, \%odds_ratio_dint, \%odds_ratio_trint ); # cria o banco de dados ou adiciona uma coluna ao banco de dados já existente

#PlotProfile($db_file, $entry_name); # Cria um arquivo PDF para a figura do perfil

#PearsonCorrelation( $db_file ); # teste estatistico utilizando Correlação de Pearson para comparar os padrões, criando uma representação da "clusterização"

print "Programa encerrado para $ARGV[0]\n";

