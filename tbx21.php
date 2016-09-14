<?php
$servername = "localhost";
$username = "root";
$password ="Pinky8109";
$database = "encode";

$conn = mysqli_connect($servername,$username,$password);
//Check Connection
if(!$conn){
    die("Connection failed: " . mysqli_connect_error());
}
//echo "connected successfully";
$sql = "USE ENCODE";
if(mysqli_query($conn,$sql)){
    //echo "Database connected";
}
else
{
    //echo "Error connecting to database";
}
        $tf = $_REQUEST["tf"];
        $limit = $_REQUEST["limit"];

        switch ($tf) {
            case "TBX21":
                # code...
                /*$filename = "data/TBX21targetgenelist.csv";
                $f = fopen($filename, "r");*/
                $sql = "SELECT * FROM TBX21TARGETGENELIST WHERE SIGNALVALUE >= $limit";
                $result = mysqli_query($conn,$sql);

                $isHeaderLine = true;
                $genelist = [];
                $list = [];
                $i = 0;
                while(($line = mysqli_fetch_assoc($result)))
                {
                    $list["gene"] = $line["V9"];
                    $list["signalvalue"] =$line["signalvalue"];
                    $genelist[$i] = $list;
                    $i = $i +1;   
                }
                

                echo json_encode($genelist);

                break;

            case "JUNB":
                # code...
                $filename = "data/JUNBtargetgenelist.csv";
                $f = fopen($filename, "r");
                $isHeaderLine = true;
                $genelist = [];
                $list = [];
                $i = 0;
                while(($line = fgets($f)) != false)
                {
                    if($isHeaderLine){
                        $isHeaderLine = false;
                    }
                    else
                    {
                    $lineArray = explode(",", $line);
                    //echo json_encode($lineArray);
                    $list["gene"] = trim($lineArray[0]);
                    $list["signalvalue"] =trim($lineArray[1]);
                    $genelist[$i] = $list;
                    $i = $i +1;
                    }
                    
                }
                fclose($f);

                echo json_encode($genelist);

                break;
            case "IRF3":
                # code...
                $filename = "data/IRF3targetgenelist.csv";
                $f = fopen($filename, "r");
                $isHeaderLine = true;
                $genelist = [];
                $list = [];
                $i = 0;
                while(($line = fgets($f)) != false)
                {
                    if($isHeaderLine){
                        $isHeaderLine = false;
                    }
                    else
                    {
                    $lineArray = explode(",", $line);
                    //echo json_encode($lineArray);
                    $list["gene"] = trim($lineArray[0]);
                    $list["signalvalue"] =trim($lineArray[1]);
                    $genelist[$i] = $list;
                    $i = $i +1;
                    }
                    
                }
                fclose($f);

                echo json_encode($genelist);

                break;
            case "SMAD5":
                # code...
                $filename = "data/SMAD5targetgenelist.csv";
                $f = fopen($filename, "r");
                $isHeaderLine = true;
                $genelist = [];
                $list = [];
                $i = 0;
                while(($line = fgets($f)) != false)
                {
                    if($isHeaderLine){
                        $isHeaderLine = false;
                    }
                    else
                    {
                    $lineArray = explode(",", $line);
                    //echo json_encode($lineArray);
                    $list["gene"] = trim($lineArray[0]);
                    $list["signalvalue"] =trim($lineArray[1]);
                    $genelist[$i] = $list;
                    $i = $i +1;
                    }
                    
                }
                fclose($f);

                echo json_encode($genelist);

                break;
            default:
                # code...
                echo "No Data Yet";
                break;
        }
    mysqli_close($conn);   
?>