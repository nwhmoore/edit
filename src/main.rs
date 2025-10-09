//https://www.cs.helsinki.fi/u/ukkonen/InfCont85.PDF

use std::fs::{File, write};
use std::io::{self, BufRead};
use std::path::Path;
use std::cmp::min;

fn main() {
    let data: Vec<String> = make_clean_fasta_data("rosalind_edit.txt");

    let ans: usize = hirschberg(&data[0],&data[1]);
    println!("Hirschberg: {}", ans);

    let ans2: Option<usize> = ukkonen(&data[0],&data[1]);
    match ans2{
        Some(dist) => println!("Ukkonen: {}", dist),
        None => println!("Ukkonen: Greater than 75% divergence")
    }
   
    //can choose ans or ans2 to write based on which algorithm you want
    write("ans.txt",ans.to_string()).expect("should write to file");
}

fn ukkonen(s1:&String, s2:&String) -> Option<usize> {
    // sort strings by length
    let a: &String;
    let b: &String;
    if s1.len() <= s2.len(){
        a = s1;
        b = s2;
    } else {
        a = s2;
        b = s1;
    }
    //a is shortest string
    //b is longer string

    let m: usize = a.len();
    let n: usize = b.len();

    //threshold! MAX DIVERGENCE
    let max_t: usize = n - (n/2)/2; //75%, can be changed in future?
    let delta: usize = 1; //minimum cost of single edit

    //initialize matrix
    let mut mat: Vec<Vec<usize>> = vec![vec![usize::MAX - 1; n+1]; m+1]; // initialize values to ~infinity
    for j in 0..n+1 {
        mat[0][j] = j;
    }
    for i in 0..m+1 {
        mat[i][0] = i;
    }

    let mut dist:Option<usize> = None;

    let mut t: usize = 1;
    let mut go: bool = true;
    while go { //while below threshold

        let mut thresh_break: bool = true; // will be used to check if all costs are above temp threshold t

        if t > max_t { // if we've increased t up to max threshold
            dist = None; // return none, past threshold
            go = false; // stop the while loop
        } else if t/delta < n-m { // if the legnths are so different they immediately violate the threshold
            t = t * 2; // increase the width of the diagonal
        } else {
            let p: usize = (t/delta - (n-m))/2; //diagonal band

            //current mat cell is [i+1][j+1]
            for i in 0..m { // for each row
                // define j loop column range here
                // restrict range to diagonal
                let col_start:usize;
                if i <= p {
                    col_start = 0;
                } else {
                    col_start = i-p;
                }
                let col_end:usize = min(n,i+(n-m)+p+1);

                thresh_break = true; // begin to check all costs vs threshold
                for j in col_start..col_end { // for each column in diagonal band
                    //evaluate d_ij from (3)
                    let cost:usize;
                    if a[i..i+1] == b[j..j+1]{ // if matching characters
                        cost = 0;
                    } else { // substitution
                        cost = delta;
                    }

                    let diag: usize = mat[i][j] + cost; // match or sub
                    let top: usize = mat[i][j+1] + delta; // indel
                    let left: usize = mat[i+1][j] + delta; // indel

                    // edit minimize cost
                    let costs: [usize; 3] = [diag,top,left];
                    let minval: Option<&usize> = costs.iter().min();
                    match minval {
                        Some(min) => mat[i+1][j+1] = *min,
                        None => println!("empty vector"),
                    }

                    if mat[i+1][j+1] <= t { // if any cost is below the temp threshold t, ok to continue
                        thresh_break = false;
                    }
                }

                if thresh_break { // break if all costs are above threshold
                    t = t * 2; // increase the temporary threshold t
                    break; // break the row for loop, start over with higher t
                }
            }
            /*
            for i in 0..m+1 {
                println!("{:?}",mat[i]); //debug
            }
            */

            if !thresh_break { // if the temp threshold t was not broken to exit the row for loop
                dist =  Some(mat[m][n]); // final element is the edit distance
                go = false; // stop while loop
            }
        }
    }

    return dist;

}

fn hirschberg(s1:&String, s2:&String) -> usize {
    //compare sizes of s1 and s2 for optimal space allocation
    let s: &String;
    let t: &String;
    /* ukkonen update broke this? maybe idk how Ordering works with str...
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
    */
    if s1.len() <= s2.len() {
        s = s1;
        t = s2;
    } else {
        s = s2;
        t = s1;
    }
    //s is shortest string
    //t is longer string

    let mut d:Vec<Vec<usize>> = vec![vec![0usize;s.len()+1];2]; // initialize columns

    // source prefixes can be transformed into empty string by
    // dropping all characters
    for i in 1..s.len()+1{
        d[0][i] = i;
    }

    // current cell is [1][i+1]
    for j in 0..t.len(){
        d[1][0] = j+1; // preallocate t prefix
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
        /* DEBUGGING
        println!("{:?}",d[0]);
        println!("{:?}",d[1]);
        println!("------");
        */
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
            let tag: &str = &line[..1]; //FASTA ID line
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