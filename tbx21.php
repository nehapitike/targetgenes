<?php
        $tf = $_REQUEST["tf"];

        switch ($tf) {
            case "TBX21":
                # code...
                $filename = "TBX21targetgenelist.csv";
                $f = fopen($filename, "r");
                $isHeaderLine = true;
                $genelist = [];
                $list = [];
                $i = 0;
                while(($line = fgets($f)) != false)
                {
                    
                    $lineArray = explode(",", $line);
                    //echo json_encode($lineArray);
                    $list["gene"] = trim($lineArray[0]);
                    $list["signalvalue"] =trim($lineArray[1]);
                    $genelist[$i] = $list;
                    $i = $i +1;
                }
                fclose($f);

                echo json_encode($genelist);

                break;
            
            default:
                # code...
                echo "No Data Yet";
                break;
        }
        
?>