// https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm

use std::fs::{File, write};
use std::io::{self, BufRead};
use std::path::Path;
use std::cmp::Ordering;

fn main() {
    let data: Vec<String> = make_clean_fasta_data("rosalind_edit.txt");

    let ans: usize = distance(&data[1],&data[0]);
    println!("{}",ans);
    write("ans.txt",ans.to_string()).expect("should write to file");
}

fn distance(s1:&String, s2:&String) -> usize {
    //compare sizes of s1 and s2 for optimal space allocation
    let s: &String;
    let t: &String;
    match s1.cmp(s2){
        Ordering::Greater => {
            s = s2;
            t = s1;
        }
        _ => {
            s = s1;
            t = s2;
        }
    }
    //s is shortest string
    //t is longer string

    let mut d:Vec<Vec<usize>> = vec![vec![0usize;s.len()+1];2]; // initialize columns

    // source prefixes can be transformed into empty string by
    // dropping all characters
    for i in 1..s.len(){
        d[0][i] = i;
    }
    // same for t prefix now
    d[1][0] = 1;

    // current cell is [1][i+1]
    for j in 0..t.len(){
        for i in 0..s.len(){

            let substitution_cost:usize;
            if s[(i)..(i+1)] == t[(j)..(j+1)]{
                substitution_cost = 0;
            } else {
                substitution_cost = 1;
            }

            let left_substitution_cost: usize = d[1][i] + 1; //deletion
            let top_substitution_cost: usize = d[0][i+1] + 1; //insertion
            let diag_substitution_cost: usize = d[0][i] + substitution_cost; //substitution

            let subs: [usize; 3] = [left_substitution_cost,top_substitution_cost,diag_substitution_cost];
            let minval: Option<&usize> = subs.iter().min();

            match minval {
                Some(min) => d[1][i+1] = *min,
                None => println!("empty vector"),
            }
        }
        d[0] = d[1].clone();

    }
    //println!("{:?}",d); //debug
    return d[1][s.len()];
}

fn make_clean_fasta_data(filepath:&str) -> Vec<String> {
    let mut data:Vec<String> = vec![String::new();0]; // final output string vector
    let mut data_line:String = String::new(); // current strand

    // File hosts.txt must exist in the current path
    if let Ok(lines) = read_lines(filepath) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines.map_while(Result::ok) {
            //println!("{}",line); //debug
            let tag = &line[..1]; //FASTA ID line
            if tag == ">" { // if ID line
                data.push(data_line); // save previous strand //MAKES EMPTY FIRST SLOT
                data_line = String::new(); // start new strand
            } else {
                data_line.push_str(&line); //push file line to current strand
            }
        }
        data.push(data_line); // push final string to strand
    } else {
        println!("Can't read file!");
    }
    data.remove(0); // remove first empty slot
    //println!("{:?}",data); //debug
    return data;
}

// The output is wrapped in a Result to allow matching on errors.
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}