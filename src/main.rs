// https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm

use std::fs::{File, write};
use std::io::{self, BufRead};
use std::path::Path;

fn main() {
    let data: Vec<String> = make_clean_fasta_data("rosalind_edit.txt");

    let ans: usize = distance(&data[0],&data[1]);
    println!("{}",ans);
    write("ans.txt",ans.to_string()).expect("should write to file");
}

// CHECK INDEXING !!!!!!!!!!!!!!!!
fn distance(s:&String, t:&String) -> usize {
    let mut d: Vec<Vec<usize>> = vec![vec![0usize;s.len()+1];t.len()+1]; // initialize matrix

    // source prefixes can be transformed into empty string by
    // dropping all characters
    for i in 1..s.len(){
        d[0][i] = i;
    }

    // same for t prefix now

    for j in 1..t.len(){
        d[j][0] = j;
    }

    // current cell is [j+1][i+1]
    for j in 0..t.len(){
        for i in 0..s.len(){

            let substitution_cost:usize;
            if s[(i)..(i+1)] == t[(j)..(j+1)]{
                substitution_cost = 0;
            } else {
                substitution_cost = 1;
            }

            let left_substitution_cost: usize = d[j+1][i] + 1; //deletion
            let top_substitution_cost: usize = d[j][i+1] + 1; //insertion
            let diag_substitution_cost: usize = d[j][i] + substitution_cost; //substitution

            let subs: [usize; 3] = [left_substitution_cost,top_substitution_cost,diag_substitution_cost];
            let minval: Option<&usize> = subs.iter().min();

            match minval {
                Some(min) => d[j+1][i+1] = *min,
                None => println!("empty vector"),
            }
        }
    }
    //println!("{:?}",d); //debug
    return d[t.len()][s.len()];
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