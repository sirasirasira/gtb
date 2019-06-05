<?php
	// args as follows
	// {aaa, abc, etc.} {1, 2, 3} {aaa, abc, etc.} {1, 2, 3}
	$cnc_arr = [-100, 0, 100];
	$from_part = strtoupper($argv[1]);
	$from_cnc = $cnc_arr[$argv[2] - 1];
	$to_part = strtoupper($argv[3]);
	$to_cnc = $cnc_arr[$argv[4] - 1];
?>
<?xml version="1.0" standalone="no"?>
<svg x="0" y="0" viewBox="150 150 300 300" version="1.1" xmlns="http://www.w3.org/2000/svg">

	<g stroke="black" stroke-width="3" fill="none" >
		<path d="M 200 200 l 100 0" stroke-dasharray="0" />
		<path d="M 300 200 l 100 0" stroke-dasharray="0" />
		<path d="M 200 400 l 100 0" stroke-dasharray="0" />
		<path d="M 300 400 l 100 0" stroke-dasharray="0" />
		<path d="M 300 300 l <?php echo $from_cnc; ?> -100" stroke-dasharray="0" />
		<path d="M 300 300 l <?php echo $to_cnc; ?> 100" stroke-dasharray="0" />
	</g>

	<g stroke="black" stroke-width="3" >
		<circle cx="200" cy="200" r="32" fill="white" />
		<text x='200' y='205' fill='black' font-size='56' stroke-width='1' text-anchor='middle' dominant-baseline='middle' ><?php echo $from_part[0] ?></text>
		<circle cx="300" cy="200" r="32" fill="white" />
		<text x='300' y='205' fill='black' font-size='56' stroke-width='1' text-anchor='middle' dominant-baseline='middle' ><?php echo $from_part[1] ?></text>
		<circle cx="400" cy="200" r="32" fill="white" />
		<text x='400' y='205' fill='black' font-size='56' stroke-width='1' text-anchor='middle' dominant-baseline='middle' ><?php echo $from_part[2] ?></text>
		<circle cx="200" cy="400" r="32" fill="white" />
		<text x='200' y='405' fill='black' font-size='56' stroke-width='1' text-anchor='middle' dominant-baseline='middle' ><?php echo $to_part[0] ?></text>
		<circle cx="300" cy="400" r="32" fill="white" />
		<text x='300' y='405' fill='black' font-size='56' stroke-width='1' text-anchor='middle' dominant-baseline='middle' ><?php echo $to_part[1] ?></text>
		<circle cx="400" cy="400" r="32" fill="white" />
		<text x='400' y='405' fill='black' font-size='56' stroke-width='1' text-anchor='middle' dominant-baseline='middle' ><?php echo $to_part[2] ?></text>
		<circle cx="300" cy="300" r="32" fill="white" />
		<text x='300' y='305' fill='black' font-size='56' stroke-width='1' text-anchor='middle' dominant-baseline='middle' >D</text>
	</g>

</svg>
